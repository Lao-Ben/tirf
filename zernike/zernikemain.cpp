#include "zernike.h"
#include <dirent.h>

using namespace std;

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

void testOnRepert(const char* s, int marge)
{
	DIR * repertoire = opendir(s);
	std::string s1(s);
	s1.append("/");

   if ( repertoire == NULL)
   {
      std::cout << "Impossible de lister le répertoire" << std::endl;
   }
   else
   {
      struct dirent * ent;

      while ( (ent = readdir(repertoire)) != NULL)
      {
      	std::string s2 = s1;
      	s2.append(ent->d_name);
      	cv::Mat src = cv::imread(s2.c_str());
			if (!src.empty())
			{
      		std::cout << ent->d_name << std::endl;
				Zernike* z = new Zernike();
				string name(ent->d_name);
				name.erase(name.size()-4, 4);
				std::vector<cv::Point2d> vect = z->getCanette(src,marge, name);
				/*if (vect.size() > 0)
				{
					for (unsigned i = 0; i < vect.size(); i++)
						std::cout << vect.at(i) << " ";
					std::cout << std::endl;
				}*/
				cv::Mat img = cv::imread(name.append("extracted.jpg"));
				double res = z->CalculateMoments(img, 5, 1);
				std::cout << res << std::endl;
			}
		}
      closedir(repertoire);
   }
}

int main(int argc, char** argv)
{
	if (atoi(argv[1]) == 1)
	{
		testOnRepert(argv[2], atoi(argv[3]));
		return 0;
	}
	cv::Mat src = cv::imread(argv[2]);
	if (src.empty())
		return -1;
		
	Zernike* z = new Zernike();
	string name(argv[2]);
	size_t found = name.find("/");
	name.erase(0,found+1);
	name.erase(name.size()-4, 4);
	std::vector<cv::Point2d> vect = z->getCanette(src,atoi(argv[3]), name);
	if (vect.size() > 0)
	{
		for (unsigned i = 0; i < vect.size(); i++)
			std::cout << vect.at(i) << " ";
		std::cout << std::endl;
	}
	cv::Mat img = cv::imread(name.append("extracted.jpg"));
	cvtColor(img, img, CV_RGB2GRAY);
	double res = z->CalculateMoments(img, 5, 1);
	std::cout << res << std::endl;	
	cv::Vec3d moy = getMoyenne(src);
	for (int i = 0; i < 3; i++)
		std::cout << moy.val[i] << " ";
	std::cout << std::endl;
	std::cout << cmplx(0, 3) << std::endl;
	return 0;
}