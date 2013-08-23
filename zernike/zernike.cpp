#include "zernike.h"
 
using namespace std; 
 
Zernike::Zernike()
{
} 

Zernike::~Zernike()
{
}

double fact(double n)
{
	if (n < 2)
		return 1;
	else
		return n*fact(n-1);
}

 double R(int p,int q, double rho)
 {
	double val=0.0;
	for(int s=0; s<=((p-abs(q))/2);s++)
	{
		double num = pow(-1,s)*fact(p-s);
		double denom = fact(s) * fact( ((p+abs(q))/2) -s ) *fact( ((p-abs(q))/2)-s );
		//double denom = fact(s) * fact( (p-2*s+q)/2) *fact( (p-2*s-q)/2);
		val += (num/denom)*pow(rho,p-2*s);
	}
	return val;
}
					
double Zernike::CalculateMoments(cv::Mat image, int n, int m)
{
    double rho, theta;
    cmplx zernike(0,0);

    for(int i=0;i<image.rows;i++)
    {
        //map y to unit circle
        for(int j=0;j<image.cols;j++)
        {
                   //map x, y  to unit ci rcle
            double y = -1 + ((2*(j+0.5))/image.cols);
            double x = -1 + ((2*(i+0.5))/image.rows);
            theta=0.0;
            rho = (double)sqrt(x*x+y*y);

            if(rho<=1.0)
            {
					if(x!=0)
						theta = atan(y/x);

					cv::Vec3b info = image.at<cv::Vec3b>(i,j);
					int val = (int) 0.1*info.val[0]+0.6*info.val[1]+0.2*info.val[2];
					cmplx pixel = val;
					cout << "i = " << i << " , j = " << j << " , x = " << x << " , y = " << y << " , rho = " << rho << " , val = " << val << " , pixel = " << pixel << endl;
					zernike+= conj( R(n,m,rho)* polar(1.0, m*theta))*pixel; //Znm = [Rnm(r)*exp(-jmO]]*f(x,y)
            }
        }
    }
    cmplx val = (n+1)/M_PI;
    zernike= zernike*val;

    double result = abs(zernike/sqrt(image.rows*image.cols)); //return normalised Zernike moment
    return result;
}