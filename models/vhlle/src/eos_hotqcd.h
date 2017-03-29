#include <vector>
class EoS ;
class Eos_hotqcd : public EoS
{
 private:
  std::vector<double> etab, stab, ptab, ttab;
	double e_min, devene;
 public:
	Eos_hotqcd(void);
	~Eos_hotqcd(void);

	virtual void eos(double e, double nb, double nq, double ns,
		double &_T, double &_mub, double &_muq, double &_mus, double &_p) ;
	virtual double p(double e, double nb, double nq, double ns) ;
        virtual void gete(double s, double& e, double nb) ;
	virtual double cs2(double e) ;
};
