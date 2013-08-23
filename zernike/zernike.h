#ifndef ZERNIKE_H
#define ZERNIKE_H

#include <complex>
#include <cmath>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>
 
 typedef std::complex<double> cmplx;

class Zernike
{
	public:
		Zernike();
		virtual ~Zernike();
		double CalculateMoments(cv::Mat image, int n, int m);
};

#endif // ZERNIKE_H