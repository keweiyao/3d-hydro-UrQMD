// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#include "event.h"

#include <iostream>
#include <algorithm>
#include <cmath>

#include <boost/program_options/variables_map.hpp>

#include <gsl/gsl_fft_complex.h>

#include "nucleus.h"

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

namespace trento {

namespace {

constexpr double TINY = 1e-12;

// Generalized mean for p > 0.
// M_p(a, b) = (1/2*(a^p + b^p))^(1/p)
inline double positive_pmean(double p, double a, double b) {
  return std::pow(.5*(std::pow(a, p) + std::pow(b, p)), 1./p);
}

// Generalized mean for p < 0.
// Same as the positive version, except prevents division by zero.
inline double negative_pmean(double p, double a, double b) {
  if (a < TINY || b < TINY)
    return 0.;
  return positive_pmean(p, a, b);
}

// Generalized mean for p == 0.
inline double geometric_mean(double a, double b) {
  return std::sqrt(a*b);
}

}  // unnamed namespace

// Determine the grid parameters like so:
//   1. Read and set step size from the configuration.
//   2. Read grid max from the config, then set the number of steps as
//      nsteps = ceil(2*max/step).
//   3. Set the actual grid max as max = nsteps*step/2.  Hence if the step size
//      does not evenly divide the config max, the actual max will be marginally
//      larger (by at most one step size).
Event::Event(const VarMap& var_map)
    : norm_(var_map["normalization"].as<double>()),
      dxy_(var_map["grid-step"].as<double>()),
      nsteps_(std::ceil(2.*var_map["grid-max"].as<double>()/dxy_)),
      xymax_(.5*nsteps_*dxy_),
      eta_steps_(var_map["eta-step"].as<int>()),
      switch_3d_(var_map["switch-3d"].as<bool>()),
      out_3d_(var_map["out-3d"].as<bool>()),
      eta_max_(var_map["eta-max"].as<double>()),
      mean_(var_map["rapidity-mean"].as<double>()),
      sigma_(var_map["rapidity-width"].as<double>()),
      skew_coeff_(var_map["rapidity-skew"].as<double>()),
      kurtosis_(var_map["rapidity-kurtosis"].as<double>()),
      K_(var_map["<pt2/mt2>"].as<double>()),
      dNdy_(boost::extents[eta_steps_]),
      TA_(boost::extents[nsteps_][nsteps_]),
      TB_(boost::extents[nsteps_][nsteps_]),
      TR_(boost::extents[nsteps_][nsteps_]),
      event_plane_list_(boost::extents[eta_steps_][2]),
      TR_3d_(boost::extents[nsteps_][nsteps_][eta_steps_]){
  // Choose which version of the generalized mean to use based on the
  // configuration.  The possibilities are defined above.  See the header for
  // more information.
  auto p = var_map["reduced-thickness"].as<double>();

  if (std::fabs(p) < TINY) {
    compute_reduced_thickness_ = [this]() {
      compute_reduced_thickness(geometric_mean);
    };
  } else if (p > 0.) {
    compute_reduced_thickness_ = [this, p]() {
      compute_reduced_thickness(
        [p](double a, double b) { return positive_pmean(p, a, b); });
    };
  } else {
    compute_reduced_thickness_ = [this, p]() {
      compute_reduced_thickness(
        [p](double a, double b) { return negative_pmean(p, a, b); });
    };
  }
}

void Event::compute(const Nucleus& nucleusA, const Nucleus& nucleusB,
                    NucleonProfile& profile) {
  // Reset npart; compute_nuclear_thickness() increments it.
  npart_ = 0;
  compute_nuclear_thickness(nucleusA, profile, TA_);
  compute_nuclear_thickness(nucleusB, profile, TB_);
  compute_reduced_thickness_();
  compute_observables();
  event_plane();
}

namespace {

// Limit a value to a range.
// Used below to constrain grid indices.
template <typename T>
inline const T& clip(const T& value, const T& min, const T& max) {
  if (value < min)
    return min;
  if (value > max)
    return max;
  return value;
}

}  // unnamed namespace

void Event::compute_nuclear_thickness(
    const Nucleus& nucleus, NucleonProfile& profile, Grid& TX) {
  // Construct the thickness grid by looping over participants and adding each
  // to a small subgrid within its radius.  Compared to the other possibility
  // (grid cells as the outer loop and participants as the inner loop), this
  // reduces the number of required distance-squared calculations by a factor of
  // ~20 (depending on the nucleon size).  The Event unit test verifies that the
  // two methods agree.

  // Wipe grid with zeros.
  std::fill(TX.origin(), TX.origin() + TX.num_elements(), 0.);

  const double r = profile.radius();

  // Deposit each participant onto the grid.
  for (const auto& nucleon : nucleus) {
    if (!nucleon.is_participant())
      continue;

    ++npart_;

    // Work in coordinates relative to (-width/2, -width/2).
    double x = nucleon.x() + xymax_;
    double y = nucleon.y() + xymax_;

    // Determine min & max indices of nucleon subgrid.
    int ixmin = clip(static_cast<int>((x-r)/dxy_), 0, nsteps_-1);
    int iymin = clip(static_cast<int>((y-r)/dxy_), 0, nsteps_-1);
    int ixmax = clip(static_cast<int>((x+r)/dxy_), 0, nsteps_-1);
    int iymax = clip(static_cast<int>((y+r)/dxy_), 0, nsteps_-1);

    // Prepare profile for new nucleon.
    profile.fluctuate();

    // Add profile to grid.
    for (auto iy = iymin; iy <= iymax; ++iy) {
      double dysq = std::pow(y - (static_cast<double>(iy)+.5)*dxy_, 2);
      for (auto ix = ixmin; ix <= ixmax; ++ix) {
        double dxsq = std::pow(x - (static_cast<double>(ix)+.5)*dxy_, 2);
        TX[iy][ix] += profile.thickness(dxsq + dysq);
      }
    }
  }
}

std::vector<double> Event::cumulent_reconstruction(double ta, double tb)
{
    double amp, arg;
    double Ym, shy, chy, y, eta; // K_ = <pt2>/<mt2>
    double k1, k2, k3, k4;
    double mu = (eta_max_ + mean_*0.5*log(ta/tb))/sigma_;
    double gamma = (ta-tb)*skew_coeff_;
    double delta = kurtosis_;
    
    std::vector<double> result_y, result_eta, eta_list;
    result_y.resize(eta_steps_);
    result_eta.resize(eta_steps_);
    eta_list.resize(eta_steps_);
    double data[2*eta_steps_];
    for(int i=0;i<eta_steps_;i++)
    {
        k1 = M_PI*(i-eta_steps_/2.0)/eta_max_*sigma_;
		k2 = k1*k1;
		k3 = k2*k1;
		k4 = k3*k1;

        amp = std::exp(-0.5*k2);
        arg = mu*k1+gamma*k3*amp;
        amp *= std::exp(-delta*k4*amp);
        
		REAL(data,i) = amp*std::cos(arg);
        IMAG(data,i) = amp*std::sin(arg);
    }
    gsl_fft_complex_radix2_forward(data, 1, eta_steps_);
    
    for(int i=0;i<eta_steps_;i++)
    {
        y = (i-eta_steps_/2.0)*2.0*eta_max_/eta_steps_;
        result_y[i] = REAL(data,i)*(2.0*static_cast<double>(i%2 == 0)-1.0);
        
        double pz = std::sinh(y);
        double pabs = std::sqrt(K_ + pz*pz);
        
        eta_list[i] = 0.5*std::log((pabs+pz)/(pabs-pz));
    }
    for(int i=0;i<eta_steps_;i++)
    {
        double eta= (i-eta_steps_/2.0)*2.0*eta_max_/eta_steps_;
	//double sh = std::sinh(eta+0.465);
	double sh = std::sinh(eta);
	double sh2 = sh*sh;
        result_eta[i] = 0.0;
        for(int ieta =0; ieta< eta_list.size()-1; ieta++)
        {
            if (eta_list[ieta] < eta && eta <= eta_list[ieta+1])
            {
                result_eta[i] = (result_y[ieta] + (result_y[ieta+1]-result_y[ieta])*(eta - eta_list[ieta])/(eta_list[ieta+1] - eta_list[ieta]))*std::sqrt(1.0+sh2)/std::sqrt(1.0/K_ + sh2);
                break;
            }
        }
    }
    return result_eta;
    
}

template <typename GenMean>
void Event::compute_reduced_thickness(GenMean gen_mean) {
  double sum = 0.;
  double ixcm = 0.;
  double iycm = 0.;

  for (int iy = 0; iy < nsteps_; ++iy) {
    for (int ix = 0; ix < nsteps_; ++ix) {
      auto t = norm_ * gen_mean(TA_[iy][ix], TB_[iy][ix]);
      TR_[iy][ix] = t;
      sum += t;
      // Center of mass grid indices.
      // No need to multiply by dxy since it would be canceled later.
      ixcm += t * static_cast<double>(ix);
      iycm += t * static_cast<double>(iy);
    }
  }

  multiplicity_ = dxy_ * dxy_ * sum;
  ixcm_ = ixcm / sum;
  iycm_ = iycm / sum;
    
//>>MODIFICATION!-------added by Weiyao-----------
  if (switch_3d_)
  {
      std::cout << "doing 3d stuff" << std::endl;
      for (int iy = 0; iy < nsteps_; ++iy) {
          for (int ix = 0; ix < nsteps_; ++ix) {
                  for (int iz = 0; iz < eta_steps_; ++iz){
                      TR_3d_[iy][ix][iz] = 0.0;
                  }
              }
          }
      for (int iz = 0; iz < eta_steps_; ++iz){
          dNdy_[iz] = 0.0;
      }
      
      std::vector<double> dSdeta;
      for (int iy = 0; iy < nsteps_; ++iy) {
          for (int ix = 0; ix < nsteps_; ++ix) {
              if(TR_[iy][ix]>0)
              {
                      dSdeta = cumulent_reconstruction(TA_[iy][ix], TB_[iy][ix]);
                      for (int iz = 0; iz < eta_steps_; ++iz){
                          if(dSdeta[iz] > 0.0 && dSdeta[int(eta_steps_/2.0)] > 0.0)
                          {
                              TR_3d_[iy][ix][iz] = TR_[iy][ix]*dSdeta[iz]/dSdeta[int(eta_steps_/2.0)];
                          }
                          else
                          {
                              TR_3d_[iy][ix][iz] = 0.0;
                          }
                          dNdy_[iz] += TR_3d_[iy][ix][iz]*dxy_*dxy_;
                      }
              }
          }
      }
  }
//<<MODIFICATION!-------added by Weiyao-----------
}

void Event::compute_observables() {
  // Compute eccentricity.

  // Simple helper class for use in the following loop.
  struct EccentricityAccumulator {
    double re = 0.;  // real part
    double im = 0.;  // imaginary part
    double wt = 0.;  // weight
    double finish() const  // compute final eccentricity
    { return std::sqrt(re*re + im*im) / std::fmax(wt, TINY); }
  } e2, e3, e4, e5;

  for (int iy = 0; iy < nsteps_; ++iy) {
    for (int ix = 0; ix < nsteps_; ++ix) {
      const auto& t = TR_[iy][ix];
      if (t < TINY)
        continue;

      // Compute (x, y) relative to the CM and cache powers of x, y, r.
      auto x = static_cast<double>(ix) - ixcm_;
      auto x2 = x*x;
      auto x3 = x2*x;
      auto x4 = x2*x2;

      auto y = static_cast<double>(iy) - iycm_;
      auto y2 = y*y;
      auto y3 = y2*y;
      auto y4 = y2*y2;

      auto r2 = x2 + y2;
      auto r = std::sqrt(r2);
      auto r4 = r2*r2;

      auto xy = x*y;
      auto x2y2 = x2*y2;

      // The eccentricity harmonics are weighted averages of r^n*exp(i*n*phi)
      // over the entropy profile (reduced thickness).  The naive way to compute
      // exp(i*n*phi) at a given (x, y) point is essentially:
      //
      //   phi = arctan2(y, x)
      //   real = cos(n*phi)
      //   imag = sin(n*phi)
      //
      // However this implementation uses three unnecessary trig functions; a
      // much faster method is to express the cos and sin directly in terms of x
      // and y.  For example, it is trivial to show (by drawing a triangle and
      // using rudimentary trig) that
      //
      //   cos(arctan2(y, x)) = x/r = x/sqrt(x^2 + y^2)
      //   sin(arctan2(y, x)) = y/r = x/sqrt(x^2 + y^2)
      //
      // This is easily generalized to cos and sin of (n*phi) by invoking the
      // multiple angle formula, e.g. sin(2x) = 2sin(x)cos(x), and hence
      //
      //   sin(2*arctan2(y, x)) = 2*sin(arctan2(y, x))*cos(arctan2(y, x))
      //                        = 2*x*y / r^2
      //
      // Which not only eliminates the trig functions, but also naturally
      // cancels the r^2 weight.  This cancellation occurs for all n.
      //
      // The Event unit test verifies that the two methods agree.
      e2.re += t * (y2 - x2);
      e2.im += t * 2.*xy;
      e2.wt += t * r2;

      e3.re += t * (y3 - 3.*y*x2);
      e3.im += t * (3.*x*y2 - x3);
      e3.wt += t * r2*r;

      e4.re += t * (x4 + y4 - 6.*x2y2);
      e4.im += t * 4.*xy*(y2 - x2);
      e4.wt += t * r4;

      e5.re += t * y*(5.*x4 - 10.*x2y2 + y4);
      e5.im += t * x*(x4 - 10.*x2y2 + 5.*y4);
      e5.wt += t * r4*r;
    }
  }

  eccentricity_[2] = e2.finish();
  eccentricity_[3] = e3.finish();
  eccentricity_[4] = e4.finish();
  eccentricity_[5] = e5.finish();
}
    
void Event::event_plane() {
    // Compute eccentricity.
        
    // Simple helper class for use in the following loop.
    struct EccentricityAccumulator {
        double re = 0.;  // real part
        double im = 0.;  // imaginary part
        double wt = 0.;  // weight
        double ex() const  // compute final eccentricity
        { return re / std::fmax(wt, TINY); }
        double ey() const  // compute event plane
        { return im / std::fmax(wt, TINY); }
    } e2;
    double ixcm_eta_, iycm_eta_, denom_;
    for(int iz = 0; iz <eta_steps_; ++iz){
        
        ixcm_eta_ = 0.0;
        iycm_eta_ = 0.0;
        denom_ = 0.0;
        for (int iy = 0; iy < nsteps_; ++iy) {
            for (int ix = 0; ix < nsteps_; ++ix) {
                const auto& t = TR_3d_[iy][ix][iz];
                ixcm_eta_ += t*static_cast<double>(ix);
                iycm_eta_ += t*static_cast<double>(iy);
                denom_ += t;
            }
        }
        ixcm_eta_/=denom_;
        iycm_eta_/=denom_;
        
        
        for (int iy = 0; iy < nsteps_; ++iy) {
            for (int ix = 0; ix < nsteps_; ++ix) {
                const auto& t = TR_3d_[iy][ix][iz];
                if (t < TINY)
                continue;
                
                // Compute (x, y) relative to the CM and cache powers of x, y, r.
                auto x = static_cast<double>(ix) - ixcm_eta_;
                auto x2 = x*x;
                
                auto y = static_cast<double>(iy) - iycm_eta_;
                auto y2 = y*y;
                
                auto r2 = x2 + y2;
                auto r = std::sqrt(r2);
                
                auto xy = x*y;
                
                // The eccentricity harmonics are weighted averages of r^n*exp(i*n*phi)
                // over the entropy profile (reduced thickness).  The naive way to compute
                // exp(i*n*phi) at a given (x, y) point is essentially:
                //
                //   phi = arctan2(y, x)
                //   real = cos(n*phi)
                //   imag = sin(n*phi)
                //
                // However this implementation uses three unnecessary trig functions; a
                // much faster method is to express the cos and sin directly in terms of x
                // and y.  For example, it is trivial to show (by drawing a triangle and
                // using rudimentary trig) that
                //
                //   cos(arctan2(y, x)) = x/r = x/sqrt(x^2 + y^2)
                //   sin(arctan2(y, x)) = y/r = x/sqrt(x^2 + y^2)
                //
                // This is easily generalized to cos and sin of (n*phi) by invoking the
                // multiple angle formula, e.g. sin(2x) = 2sin(x)cos(x), and hence
                //
                //   sin(2*arctan2(y, x)) = 2*sin(arctan2(y, x))*cos(arctan2(y, x))
                //                        = 2*x*y / r^2
                //
                // Which not only eliminates the trig functions, but also naturally
                // cancels the r^2 weight.  This cancellation occurs for all n.
                //
                // The Event unit test verifies that the two methods agree.
                e2.re += t * (y2 - x2);
                e2.im += t * 2.*xy;
                e2.wt += t * r2;
            }
        }
        event_plane_list_[iz][0] = e2.ex();
        event_plane_list_[iz][1] = e2.ey();
        
    }
}
}  // namespace trento
