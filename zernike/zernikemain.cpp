#include "zernike.h"

int main(int argc, char** argv)
{
	cv::Mat src = cv::imread(argv[1]);
	if (src.empty())
		return -1;
		
	Zernike* z = new Zernike();
	double res = z->CalculateMoments(src, 5, 1);	
	std::cout << res << std::endl;
	return 0;
}