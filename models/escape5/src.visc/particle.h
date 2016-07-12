class ParticlePDG2 ;

class Particle {
public:
  double px, py, pz, e ;
  double x, y, z, t ;
  ParticlePDG2 *def ;
  int mid ;
  Particle(double X, double Y, double Z, double T,
    double Px, double Py, double Pz, double E, ParticlePDG2* def1, int Mid):
    px(Px), py(Py), pz(Pz), e(E), x(X), y(Y), z(Z), t(T), def(def1), mid(Mid) {} ;
  ~Particle() {} ;
} ;
