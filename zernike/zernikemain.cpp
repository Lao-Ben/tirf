#include "zernike.h"

cv::Vec3d getMoyenne(cv::Mat image)
{
	int rows = image.rows;
	int cols = image.cols;
	
	cv::Vec3d res(0,0,0);	
	int tot = 0;
	
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
		{
			cv::Vec3d info = image.at<cv::Vec3b>(i,j);
			if ((info.val[0] != 0 || info.val[1] != 0 || info.val[2] != 0) && (info.val[0] != 255 || info.val[1] != 255 || info.val[2] != 255))
			{				
				res = res + info;
				tot++;
			}
		}
	cv::Vec3d ret(res.val[0]/tot,res.val[1]/tot,res.val[2]/tot);
	return ret;
}

int main(int argc, char** argv)
{
	cv::Mat src = cv::imread(argv[1]);
	if (src.empty())
		return -1;
		
	Zernike* z = new Zernike();
	std::vector<cv::Point2d> vect = z->getCanette(src);
	/*if (vect.size() > 0)
	{
		for (unsigned i = 0; i < vect.size(); i++)
			std::cout << vect.at(i) << " ";
		std::cout << std::endl;
	}*/
	cv::Mat img = cv::imread("extracted.jpg");
	double res = z->CalculateMoments(img, 5, 1);
	std::cout << res << std::endl;	
	cv::Vec3d moy = getMoyenne(src);
	for (int i = 0; i < 3; i++)
		std::cout << moy.val[i] << " ";
	std::cout << std::endl;
	return 0;
}