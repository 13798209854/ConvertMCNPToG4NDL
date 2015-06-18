#include "../include/AngularEnergyDist.hh"

AngularEnergyDist::AngularEnergyDist()
{
    //ctor
}

AngularEnergyDist::~AngularEnergyDist()
{
    //dtor
}

double AngularEnergyDist::Interpolate(int aScheme, double x, double x1, double x2, double y1, double y2) const
{
  double result(0);
  int theScheme = aScheme;
  theScheme = theScheme%7;
  switch(theScheme)
  {
    case 1:
      //080809
      //result = Histogram(x, x1, x2, y1, y2);
      result = LinearLinear(x, x1, x2, y1, y2);
      break;
    case 2:
      result = LinearLinear(x, x1, x2, y1, y2);
      break;
    case 3:
      result = LinearLogarithmic(x, x1, x2, y1, y2);
      break;
    case 4:
      result = LogarithmicLinear(x, x1, x2, y1, y2);
      break;
    case 5:
      result = LogarithmicLogarithmic(x, x1, x2, y1, y2);
      break;
    default:
      cout << "Error: Unrecognized scheme = "<<theScheme<<endl;
      break;
  }
  return result;
}

inline double AngularEnergyDist::
Histogram(double , double , double , double y1, double ) const
{
  double result;
  result = y1;
  return result;
}

inline double AngularEnergyDist::
LinearLinear(double x, double x1, double x2, double y1, double y2) const
{
  double slope=0, off=0;
  if(x2-x1==0) return (y2+y1)/2.;
  slope = (y2-y1)/(x2-x1);
  off = y2-x2*slope;
  double y = x*slope+off;
  return y;
}

inline double AngularEnergyDist::
LinearLogarithmic(double x, double x1, double x2, double y1, double y2) const
{
  double result;
  if(x==0) result = y1+y2/2.;
  else if(x1==0) result = y1;
  else if(x2==0) result = y2;
  else result = LinearLinear(std::log(x), std::log(x1), std::log(x2), y1, y2);
  return result;
}

inline double AngularEnergyDist::
LogarithmicLinear(double x, double x1, double x2, double y1, double y2) const
{
  double result;
  if(y1==0||y2==0) result = 0;
  else
  {
    result = LinearLinear(x, x1, x2, std::log(y1), std::log(y2));
    result = std::exp(result);
  }
  return result;
}

inline double AngularEnergyDist::
LogarithmicLogarithmic(double x, double x1, double x2, double y1, double y2) const
{
  double result;
  if(x==0) result = y1+y2/2.;
  else if(x1==0) result = y1;
  else if(x2==0) result = y2;
  if(y1==0||y2==0) result = 0;
  else
  {
    result = LinearLinear(std::log(x), std::log(x1), std::log(x2), std::log(y1), std::log(y2));
    result = std::exp(result);
  }
  return result;
}
