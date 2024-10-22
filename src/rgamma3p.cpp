#include <RcppEigen.h>

#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE
#define pi 3.141592653589793
#define s2pi 2.506628274631

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppEigen)]]

#include "factorial.h"
#include "pbdv.h"


template <typename T> int sign (T val) {
  return (T(0) < val) - (val < T(0));
}

// [[Rcpp::export]]

double rg3p_c1 (const double& a, const double& b) 
{
  double r = std::abs(b) / a;
  int rs;
  
  if (r <= 0.2) {
    rs = 1;
  } else if (r > 0.2 && r <= 0.4) {
    rs = 2;
  } else if (r > 0.4 && r <= 0.6) {
    rs = 3;
  } else if (r > 0.6 && r <= 0.8) {
    rs = 4;
  } else if (r > 0.8 && r <= 1.0) {
    rs = 5;
  } else if (r > 1.0 && r <= 1.5) {
    rs = 6;
  } else if (r > 1.5 && r <= 2.0) {
    rs = 7;
  } else if (r > 2.0 && r <= 2.5) {
    rs = 8;
  } else if (r > 2.5 && r <= 3.0) {
    rs = 9;
  } else if (r > 3.0 && r <= 3.5) {
    rs = 10;
  } else if (r > 3.5 && r <= 4.0) {
    rs = 11;
  } else if (r > 4.0 && r <= 4.5) {
    rs = 12;
  } else if (r > 4.5 && r <= 5.0) {
    rs = 13;
  } else {
    rs = 0;
  }
  
  double out = 0.0;
  
  if (b < 0.0) 
  {
    // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // case b < 0
    
    if (r >= 5.0) 
    {
      // Approximation to Gamma distribution
      if (r >= 20.0) 
      {
        // approximated parameters
        out = R::rgamma(2.0, 1.0 / std::abs(b));
        return out;
      } 
      else 
      {
        // exact parameters
        const double ts = std::sqrt(2.0) * a;
        const double t = -b / ts;
        
        const double d21 = pbdv(-3.0, t) / pbdv(-2.0, t);
        const double d31 = (1.0 - t * d21) / 3.0;
        
        const double mu = 2.0 / ts * d21;
        const double s2 = 6.0 / (ts * ts) * d31 - mu * mu;  
        
        const double u = mu / s2;
        out = R::rgamma(u * mu, 1.0 / u);
        return out;
      }
    } 
    else 
    {
      // Exact sampling
      double xs1 = 0.0, xs2 = 0.0, ds = 0.0, bs = 0.0, cs = 0.0;

      switch(rs) 
      {
      case 1:
        xs1 = -1.5305156071382; xs2 = 0.912689119744102;
        ds = 0.655348065027528; bs = -0.4964561; cs = 0.690138935200033;
        break; 
      case 2:
        xs1 = -1.52631273299577; xs2 = 0.929115435735049;
        ds = 0.659234853596739; bs = -0.5028346; cs = 0.73468663501275;
        break; 
      case 3:
        xs1 = -1.52162340688057; xs2 = 0.94628704376105;
        ds = 0.662997556152702; bs = -0.5085152; cs = 0.779040199076203;
        break; 
      case 4:
        xs1 = -1.51659405713692; xs2 = 0.964115091624826;
        ds = 0.666656721660187; bs = -0.5135407; cs = 0.822936278311114;
        break; 
      case 5:
        xs1 = -1.51135018545505; xs2 = 0.982514613033376;
        ds = 0.670231423536598; bs = -0.5179581; cs = 0.86613962768219;
        break; 
      case 6:
        xs1 = -1.49794734172668; xs2 = 1.00140519421345;
        ds = 0.678903806481405; bs = -0.5240614; cs = 0.969856251460308;
        break; 
      case 7:
        xs1 = -1.48498179745043; xs2 = 1.05029781972212;
        ds = 0.687345894831465; bs = -0.530846; cs = 1.06573016608705;
        break; 
      case 8:
        xs1 = -1.47308147718974; xs2 = 1.10076692992126;
        ds = 0.695664093100522; bs = -0.5353669; cs = 1.15257389067033;
        break; 
      case 9:
        xs1 = -1.46251252652624; xs2 = 1.15197785672063;
        ds = 0.703896084751363; bs = -0.5382237; cs = 1.23008749716618;
        break; 
      case 10:
        xs1 = -1.45332293238867; xs2 = 1.20327937182714;
        ds = 0.712036465716239; bs = -0.5398922; cs = 1.29857308211284;
        break; 
      case 11:
        xs1 = -1.44544184030771; xs2 = 1.25417822831087;
        ds = 0.720058127965459; bs = -0.5407328; cs = 1.35868596210447;
        break; 
      case 12:
        xs1 = -1.43874240881106; xs2 = 1.30431194740893;
        ds = 0.727926692044799; bs = -0.5410095; cs = 1.41124992334575;
        break; 
      case 13:
        xs1 = -1.43307882982896; xs2 = 1.35342319244879;
        ds = 0.735608926519893; bs = -0.5409117; cs = 1.45713471624584;
        break;
      }
      
      const double ts = std::sqrt(2.0) * a;
      const double t = -b / ts;
      const double d1 = pbdv(-2, t);
      const double d21 = pbdv(-3, t) / d1;
      const double d31 = (1.0 - t * d21);
      
      const double mu = 2.0 * d21 / ts;
      const double s2 = 2.0 * d31 / (ts * ts) - mu * mu; 
      
      const double s = std::sqrt(s2);
      const double w = 1.0 / (ts * s);
      
      double u = R::rnorm(0.0, w);
      double kf = 0.0, kg = 0.0, q = 0.0;
      
      out = u * s + mu;
      
      if ((xs1 <= u) && (u <= xs2)) 
      {  
        // :::::::::::: 
        // STEP 1
        
        return out;
      } 
      else 
      {
        // :::::::::::: 
        // STEP 2
        
        kf = s * ts * ts / (std::exp(0.25 * t * t) * d1);
        kg = w * s2pi;
        
        if (u > -mu/s) 
        {
          q = kf * kg * out * std::exp(out * (b - ts * ts * mu) + a * a * mu * mu);
          if (R::runif(0, 1) <= q) 
            return out;
        }
      }
      
      // STEP 3
      int i = 0;
      double e, x, f, g, dlap;
      while (i < 1e+16) 
      {
        i += 1;
        e = R::rexp(1);
        u = R::runif(0, 1);
        u = u + u - 1.0;
        
        if (u < 0.0) 
        {
          x = bs - ds * e;
        } 
        else 
        {
          x = bs + ds * e;
        }
        out = x * s + mu;
        
        if ((xs1 <= x) && (x <= xs2)) 
        {
          f = kf * out * std::exp(-a * a * out * out + b * out);
          g = std::exp(-0.5 * x * x / (w * w)) / kg;
          dlap = cs / s2pi * std::exp(-std::abs(x - bs) / ds);
          
          if (std::abs(u) <= ((f - g) / dlap)) 
            return out;
        }
      }
    }
  } 
  else 
  {
    // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // case b > 0
    
    if (r >= 5.0) 
    {
      // Approximation to Normal distribution
      if (r >= 30.0)
      {
        // approximated parameters
        const double q = 2.0 * a * a;
        const double mu = b / q;
        const double s2 = 1.0 / q;

        out = R::rnorm(mu, sqrt(s2));
        while (out <= 0.0)
          out = R::rnorm(mu, sqrt(s2));
        
        return out;
      }
      else
      {
        // exact parameters
        const double ts = std::sqrt(2.0) * a;
        const double t = -b / ts;

        const double d21 = pbdv(-3.0, t) / pbdv(-2.0, t);
        const double d31 = (1.0 - t * d21) / 3.0;

        const double mu = 2.0 / ts * d21;
        const double s2 = 6.0 / (ts * ts) * d31 - mu * mu;

        out = R::rnorm(mu, sqrt(s2));
        while (out <= 0.0)
          out = R::rnorm(mu, sqrt(s2));

        return out;
      }
    } 
    else 
    {
      // Exact sampling
      double xs1 = 0.0, xs2 = 0.0, ds = 0.0, bs = 0.0, cs = 0.0;
      
      switch(rs) 
      {
      case 1:
        xs1 = -1.53406499365263; xs2 = 0.897099868956042;
        ds = 0.651316106787244; bs = -0.4893424; cs = 0.645683883557984;
        break; 
      case 2:
        xs1 = -1.53677503595901; xs2 = 0.882441053489848;
        ds = 0.647118639185317; bs = -0.481464; cs = 0.601625204378373;
        break; 
      case 3:
        xs1 = -1.53844474387079; xs2 = 0.868805896139059;
        ds = 0.642737657643309; bs = -0.4728009; cs = 0.558275397563083;
        break; 
      case 4:
        xs1 = -1.5378220009327; xs2 = 0.856285253741137;
        ds = 0.638159318660735; bs = -0.4633445; cs = 0.515946250181509;
        break; 
      case 5:
        xs1 = -1.53511863220306; xs2 = 0.84496488884044;
        ds = 0.633375865516132; bs = -0.4530994; cs = 0.474938553341791;
        break; 
      case 6:
        xs1 = -1.52003476398805; xs2 = 0.822392749748821;
        ds = 0.628387467204558; bs = -0.432877; cs = 0.435531694260361;
        break; 
      case 7:
        xs1 = -1.49217848955996; xs2 = 0.808571949747569;
        ds = 0.615108869353087; bs = -0.4011468; cs = 0.345547749251269;
        break; 
      case 8:
        xs1 = -1.45284887462643; xs2 = 0.803454556472932;
        ds = 0.601074984418797; bs = -0.3663053; cs = 0.26968974729252;
        break; 
      case 9:
        xs1 = -1.4061695417862; xs2 = 0.803454556472932;
        ds = 0.587071622998585; bs = -0.3300687; cs = 0.208720977879743;
        break; 
      case 10:
        xs1 = -1.35770976619668; xs2 = 0.80591579713146;
        ds = 0.574019376896896; bs = -0.2944748; cs = 0.161643333219136;
        break; 
      case 11:
        xs1 = -1.31238194093465; xs2 = 0.813896651189013;
        ds = 0.562675775766095; bs = -0.2614344; cs = 0.12630407612564;
        break; 
      case 12:
        xs1 = -1.27303663426414; xs2 = 0.824969842044926;
        ds = 0.553414886421639; bs = -0.2322547; cs = 0.100143330261176;
        break; 
      case 13:
        xs1 = -1.24039674144834; xs2 = 0.83703314448015;
        ds = 0.546199377716538; bs = -0.2074166; cs = 0.080779310328355;
        break;
      }
      
      const double ts = std::sqrt(2.0) * a;
      const double t = -b / ts;
      const double d1 = pbdv(-2, t);
      const double d21 = pbdv(-3, t) / d1;
      const double d31 = (1.0 - t * d21);
      
      const double mu = 2.0 * d21 / ts;
      const double s2 = 2.0 * d31 / (ts * ts) - mu * mu; 
      
      const double s = std::sqrt(s2);
      const double w = 1.0 / (ts * s);
      
      double u = R::rnorm(0.0, w);
      double kf = 0.0, kg = 0.0, q = 0.0;
      
      out = u * s + mu;
      
      if ((xs1 <= u) && (u <= xs2)) 
      {  
        // :::::::::::: 
        // STEP 1
        
        return out;
      } 
      else 
      {
        // :::::::::::: 
        // STEP 2
        
        kf = s * ts * ts / (std::exp(0.25 * t * t) * d1);
        kg = w * s2pi;
        
        if (u > -mu/s) 
        {
          q = kf * kg * out * std::exp(out * (b - ts * ts * mu) + a * a * mu * mu);
          if (R::runif(0, 1) <= q) 
            return out;
        }
      }
      
      // STEP 3
      int i = 0;
      double e, x, f, g, dlap;
      while (i < 1e+16) 
      {
        i += 1;
        e = R::rexp(1);
        u = R::runif(0, 1);
        u = u + u - 1.0;
        
        if (u < 0.0) 
        {
          x = bs - ds * e;
        } 
        else 
        {
          x = bs + ds * e;
        }
        out = x * s + mu;
        
        if ((xs1 <= x) && (x <= xs2)) 
        {
          f = kf * out * std::exp(-a * a * out * out + b * out);
          g = std::exp(-0.5 * x * x / (w * w)) / kg;
          dlap = cs / s2pi * std::exp(-std::abs(x - bs) / ds);
          
          if (std::abs(u) <= ((f - g) / dlap)) 
            return out;
        }
      }
    }
  }
  
  Rcpp::Rcout << "\n Warning rg3p_c1: returned value 0 \n";
  return 0.0;
}


// [[Rcpp::export]]

double rg3p_approx (const double& a, const double& b, const int& c) 
{
  if ((c > 14) || (b > 0.0)) 
  {
    const double s2a  = std::sqrt(2.0) * a;
    const double bs2a = -b / s2a;
    const double bs2a2 = bs2a * bs2a;
    
    const double q = -bs2a + 0.5 * (bs2a + std::sqrt(bs2a2 + 4 * c + 2));
    const double w = -bs2a + 0.5 * (bs2a + std::sqrt(bs2a2 + 4 * c + 6));
    
    const double mu = q / s2a;
    const double s2 = q * (w - q) / (s2a * s2a);
    
    double out = R::rnorm(mu, std::sqrt(s2));
    while (out <= 0.0) 
      out = R::rnorm(mu, std::sqrt(s2));  

    return out;
  }
  else
  {
    double r = std::abs(b) / a;
    
    if (r >= 15.0) {
      // approximated parameters
      const double u = std::abs(b);
      const double q = c + 1;
      
      const double out = R::rgamma(q, 1 / u);
      return out;
      
    } 
    else 
    {
      // exact parameters
      const double ts = std::sqrt(2.0) * a;
      const double t = -b / ts;
      
      const double d21 = pbdv(-c-2, t) / pbdv(-c-1, t);
      const double d31 = (1.0 - t * d21) / (c+2);
      
      const double mu = (c+1) * d21 / ts;
      const double s2 = (c+1) * (c+2) / (ts * ts) * d31 - mu * mu;
      // const double mu = (c+1) * pbdv(-c-2, t) / (d1 * ts);
      // const double s2 = (c+1) * (c+2) / (ts*ts) * pbdv(-c-3, t) / d1 - mu * mu;  
      
      const double u = mu / s2;
      const double q = u * mu;
      
      const double out = R::rgamma(q, 1 / u);
      return out;
    } 
  }
  
  Rcpp::Rcout << "\n Warning rg3p_approx: returned value 0 \n";
  return 0.0;
}


// end file