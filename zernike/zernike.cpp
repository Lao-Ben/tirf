#include "zernike.h"
 
using namespace std; 
using namespace cv;
 
int nbPointAround(Mat image, Point2d pt, cv::Vec3b col, int marge); 
Mat detectObjectWithColor(Mat image, Vec3b scalColor, int margeR, int margeG, int margeB, Mat img);
 
Zernike::Zernike()
{
} 

Zernike::~Zernike()
{
}

bool goodColor(cv::Vec3b pix, cv::Vec3b col, int margeR, int margeG, int margeB)
{
	if(pix.val[0] > (col.val[0] - margeB) && (pix.val[0] < (col.val[0] + margeB)))
		if(pix.val[1] > (col.val[1] - margeG) && (pix.val[1] < (col.val[1] + margeG)))
			if(pix.val[2] > (col.val[2] - margeR) && (pix.val[2] < (col.val[2] + margeR)))
				return true;
	return false;
}

bool goodColor(cv::Vec3b pix, cv::Vec3b col, int marge)
{
	int margeR = marge;
   int margeG = marge;
   int margeB = marge;
	if(pix.val[0] > (col.val[0] - margeB) && (pix.val[0] < (col.val[0] + margeB)))
		if(pix.val[1] > (col.val[1] - margeG) && (pix.val[1] < (col.val[1] + margeG)))
			if(pix.val[2] > (col.val[2] - margeR) && (pix.val[2] < (col.val[2] + margeR)))
				return true;
	return false;
}

vector<double> getFinalVect(vector<double> vect, bool min)
{
	vector<double> res;
	if (min)
	{
		res.push_back(vect.at(0));
		for (unsigned int i=1; i < vect.size(); i++)
			if (vect.at(i)-vect.at(i-1) > 100)
				res.push_back(vect.at(i));
	}
	else
	{
		res.push_back(vect.at(vect.size()-1));
		for (unsigned int i=vect.size()-2; i > 0; i--)
			if (abs(vect.at(i)-vect.at(i-1)) > 100)
				res.push_back(vect.at(i));
		reverse(res.begin(),res.end());
	}
	return res;
}

int nbPointAround(Mat image, Point2d pt, cv::Vec3b col, int marge)
{
	double x = pt.x;
	double y = pt.y;
	
	int res = 0;
	for(unsigned int i=x-6;i <= x+6; i++)
	{
		for (unsigned int j=y-6; j <= y+6; j++)
		{
			if (i >= 0 && i < image.rows && j >= 0 && j < image.cols)
				if (goodColor(image.at<Vec3b>(i,j), col, marge))
					res = res + 1;
		}
	}
	return res;
}

int nbPointAround(Mat image, Point2d pt, cv::Vec3b col, int margeR, int margeG, int margeB)
{
	double x = pt.x;
	double y = pt.y;
	
	int res = 0;
	for(unsigned int i=x-6;i <= x+6; i++)
	{
		for (unsigned int j=y-6; j <= y+6; j++)
		{
			if (i >= 0 && i < image.rows && j >= 0 && j < image.cols)
				if (goodColor(image.at<Vec3b>(i,j), col, margeR, margeG, margeB))
					res = res + 1;
		}
	}
	return res;
}

vector<vector<Point2d> > getRects(Mat image, vector<Point2d> vect, Vec3b scalColor, int marge)
{
	vector<vector<Point2d> > res;
	vector<double> vectMinX;
	vector<double> vectMinY;
	vector<double> vectMaxX;
	vector<double> vectMaxY;
	for (unsigned int i=0; i < vect.size(); i++)
	{
		Point2d pt = vect.at(i);
		//if(pt.y == 0 || vectMinX.size() == 0 || pt.y < vectMinX.at(vectMinX.size()-1))
			vectMinX.push_back(pt.x);
		//if(pt.y == image.rows || vectMaxX.size() == 0 || pt.y > vectMaxX.at(vectMaxX.size()-1))
			vectMaxX.push_back(pt.x);
		//if(pt.x == 0 || vectMinY.size() == 0 || pt.x < vectMinY.at(vectMinY.size()-1))
			vectMinY.push_back(pt.y);
		//if(pt.x == image.cols || vectMaxY.size() == 0 || pt.x > vectMaxY.at(vectMaxY.size()-1))
			vectMaxY.push_back(pt.y);
	}
	std::vector<double>::iterator it;
	sort(vectMinX.begin(),vectMinX.end());
  	it = unique(vectMinX.begin(),vectMinX.end());
  	vectMinX.resize(distance(vectMinX.begin(),it));
	sort(vectMinY.begin(),vectMinY.end());
	it = unique(vectMinY.begin(),vectMinY.end());
  	vectMinY.resize(distance(vectMinY.begin(),it));
	sort(vectMaxX.begin(),vectMaxX.end());
	it = unique(vectMaxX.begin(),vectMaxX.end());
  	vectMaxX.resize(distance(vectMaxX.begin(),it));
	sort(vectMaxY.begin(),vectMaxY.end());
	it = unique(vectMaxY.begin(),vectMaxY.end());
  	vectMaxY.resize(distance(vectMaxY.begin(),it));
  	vector<double> vectMinXfinal = getFinalVect(vectMinX, true);
	vector<double> vectMinYfinal = getFinalVect(vectMinY, true);
	vector<double> vectMaxXfinal = getFinalVect(vectMaxX, false);
	vector<double> vectMaxYfinal = getFinalVect(vectMaxY, false);
	for (unsigned int i=0; i < vectMinXfinal.size(); i++)
	{
		vector<Point2d> v;
		v.push_back(Point2d(vectMinYfinal.at(i), vectMinXfinal.at(i)));
		v.push_back(Point2d(vectMinYfinal.at(i),vectMaxXfinal.at(i)));
		v.push_back(Point2d(vectMaxYfinal.at(i),vectMinXfinal.at(i)));
		v.push_back(Point2d(vectMaxYfinal.at(i),vectMaxXfinal.at(i)));
		for (unsigned j = 0; j < v.size(); j++)
			std::cout << v.at(j) << " ";
		std::cout << std::endl;
		res.push_back(v);
	}
	return res;
}

std::vector<cv::Point2d> getRect(Mat image, std::vector<cv::Point2d> vect, cv::Vec3b col, int margeR, int margeG, int margeB)
{
	vector<Point2d> res;
	double minX = -1;
	double minY = -1;
	double maxX = -1;
	double maxY = -1;
	for (unsigned int i=0; i < vect.size(); i++)
	{
		Point2d pt = vect.at(i);
		int nbArround = nbPointAround(image, pt, col, margeR, margeG, margeB);
		if (nbArround >= 165)
		{
			if (minX == -1 || minX > pt.x)
				minX = pt.x;
			if (minY == -1 || minY > pt.y)
				minY = pt.y;
			if (maxX == -1 || maxX < pt.x)
				maxX = pt.x;
			if (maxY == -1 || maxY < pt.y)
				maxY = pt.y;
		}
	}
	res.push_back(Point2d(minY,minX));
	res.push_back(Point2d(minY,maxX));
	res.push_back(Point2d(maxY,minX));
	res.push_back(Point2d(maxY,maxX));
	return res;
}

std::vector<cv::Point2d> getRect(Mat image, std::vector<cv::Point2d> vect, cv::Vec3b col, int marge)
{
	vector<Point2d> res;
	double minX = -1;
	double minY = -1;
	double maxX = -1;
	double maxY = -1;
	for (unsigned int i=0; i < vect.size(); i++)
	{
		Point2d pt = vect.at(i);
		int nbArround = nbPointAround(image, pt, col, marge);
		if (nbArround >= 165)
		{
			if (minX == -1 || minX > pt.x)
				minX = pt.x;
			if (minY == -1 || minY > pt.y)
				minY = pt.y;
			if (maxX == -1 || maxX < pt.x)
				maxX = pt.x;
			if (maxY == -1 || maxY < pt.y)
				maxY = pt.y;
		}
	}
	res.push_back(Point2d(minY,minX));
	res.push_back(Point2d(minY,maxX));
	res.push_back(Point2d(maxY,minX));
	res.push_back(Point2d(maxY,maxX));
	return res;
}

void rotate(cv::Mat& src, double angle, cv::Mat& dst)
{
    int len = std::max(src.cols, src.rows);
    cv::Point2f pt(len/2., len/2.);
    cv::Mat r = cv::getRotationMatrix2D(pt, angle, 1.0);

    cv::warpAffine(src, dst, r, cv::Size(dst.cols, dst.rows));
}

double percentGoodColor(Mat image, Vec3b col, int marge)
{
	int rows = image.rows;
	int cols = image.cols;
	double res = 0;
	
	for (unsigned int i = 0; i < rows; i++)
		for (unsigned int j = 0; j < cols; j++)
		{
			Vec3b pix = image.at<Vec3b>(i,j);
			if(goodColor(pix, col, marge))
				res = res + 1;
		}
	res = 100 * res / (rows*cols);
	return res;
}

double percentGoodColor(Mat image, Vec3b col, int margeR, int margeG, int margeB)
{
	int rows = image.rows;
	int cols = image.cols;
	double res = 0;
	
	for (unsigned int i = 0; i < rows; i++)
		for (unsigned int j = 0; j < cols; j++)
		{
			Vec3b pix = image.at<Vec3b>(i,j);
			if(goodColor(pix, col, margeR, margeG, margeB))
				res = res + 1;
		}
	res = 100 * res / (rows*cols);
	return res;
}

Mat extractRect(Mat image, std::vector<cv::Point2d> rect, Vec3b scalColor, int margeR, int margeG, int margeB, Mat img)
{
	double minX = rect.at(0).x;
	double minY = rect.at(0).y;
	double maxX = rect.at(3).x;
	double maxY = rect.at(3).y; 
	Mat M(maxY-minY+1,maxX-minX+1, image.type());
	for (unsigned int i = minX; i <= maxX;i++)
		for (unsigned int j = minY; j <= maxY;j++)
			M.at<Vec3b>(j-minY,i-minX) = image.at<Vec3b>(j,i);
	 if (M.rows > M.cols)
    {
    	Mat tmp(maxX-minX+1,maxY-minY+1, M.type());
    	rotate(M,-90, tmp);
    	M = tmp;
    }
    double percent = percentGoodColor(M, scalColor, margeR, margeG, margeB);
   if (percent <= 15 && percent > 0)
   	M = detectObjectWithColor(image, scalColor, margeR, margeG-2, margeB-2, img);
   else
   {
    	line(img, rect.at(0), rect.at(1), Scalar(0,255,0), 1, 8, 0);
    	line(img, rect.at(0), rect.at(2), Scalar(0,255,0), 1, 8, 0);
    	line(img, rect.at(1), rect.at(3), Scalar(0,255,0), 1, 8, 0);
    	line(img, rect.at(2), rect.at(3), Scalar(0,255,0), 1, 8, 0);
   }
   return M;
}

void extractRect(Mat image, std::vector<cv::Point2d> rect, string name, Vec3b scalColor, int marge)
{
	double minX = rect.at(0).x;
	double minY = rect.at(0).y;
	double maxX = rect.at(3).x;
	double maxY = rect.at(3).y; 
	Mat M(maxY-minY+1,maxX-minX+1, image.type());
	for (unsigned int i = minX; i <= maxX;i++)
		for (unsigned int j = minY; j <= maxY;j++)
			M.at<Vec3b>(j-minY,i-minX) = image.at<Vec3b>(j,i);
	 if (M.rows > M.cols)
    {
    	Mat tmp(maxX-minX+1,maxY-minY+1, M.type());
    	rotate(M,-90, tmp);
    	M = tmp;
    }
    double percent = percentGoodColor(M, scalColor, marge);
   if (percent <= 15 && percent > 0)
   	M = detectObjectWithColor(image, scalColor, marge, marge-2, marge-2, image);
   else
   {
    	line(image, rect.at(0), rect.at(1), Scalar(0,255,0), 1, 8, 0);
    	line(image, rect.at(0), rect.at(2), Scalar(0,255,0), 1, 8, 0);
    	line(image, rect.at(1), rect.at(3), Scalar(0,255,0), 1, 8, 0);
    	line(image, rect.at(2), rect.at(3), Scalar(0,255,0), 1, 8, 0);
   }
   Mat img(140,280, M.type());
   resize(M,img,img.size(),0,0,INTER_CUBIC);
	imwrite(name.append("extracted.jpg"), img);
}

Mat detectObjectWithColor(Mat image, Vec3b scalColor, int margeR, int margeG, int margeB, Mat img)
{
    Vec3b scalPix;
    Vec3b scalPainter;
    std::vector<cv::Point2d> points;
    
    Mat img_object = image.clone();
 
    for(int x=0; x<img_object.rows; x++)
    {
    	for(int y=0; y<img_object.cols; y++)
		{
            //On récupère le pixel p(x;y)
            	Vec3b scalPix = img_object.at<Vec3b>(x,y);
 
                //On teste sa couleur
                if (goodColor(scalPix, scalColor, margeR))
                {
                	scalPainter = scalPix;
                	points.push_back(Point2d(x,y));
                }
                else
                {
                    scalPainter.val[0] = 0;
                    scalPainter.val[1] = 0;
                    scalPainter.val[2] = 0;
                }
         	img_object.at<Vec3b>(x,y) = scalPainter;
        }
    }
    vector<Point2d> rect = getRect(image, points, scalColor, margeR, margeG, margeB);
    Mat mat = extractRect(image, rect, scalColor, margeR, margeG, margeB, img);
    return mat;
}

vector<Point2d> detectObjectWithColor(Mat image, Vec3b scalColor, int marge, string name)
{
    int margeR = marge;
    int margeG = marge;
    int margeB = marge;
 
    Vec3b scalPix;
    Vec3b scalPainter;
    std::vector<cv::Point2d> points;
    
    Mat img_object = image.clone();
 
    for(int x=0; x<img_object.rows; x++)
    {
    	for(int y=0; y<img_object.cols; y++)
		{
            //On récupère le pixel p(x;y)
            	Vec3b scalPix = img_object.at<Vec3b>(x,y);
 
                //On teste sa couleur
                if (goodColor(scalPix, scalColor, marge))
                {
                	scalPainter = scalPix;
                	points.push_back(Point2d(x,y));
                }
                else
                {
                    scalPainter.val[0] = 0;
                    scalPainter.val[1] = 0;
                    scalPainter.val[2] = 0;
                }
         	img_object.at<Vec3b>(x,y) = scalPainter;
        }
    }
    vector<Point2d> rect = getRect(image, points, scalColor, marge);
    extractRect(image, rect, name, scalColor, marge);
    /*vector<vector<Point2d> > rects = getRects(image, points, scalColor, marge);
    for(unsigned int i = 0; i < rects.size(); i++)
    {
    	line(image, rects.at(i).at(0), rects.at(i).at(1), Scalar(0,255,0), 1, 8, 0);
    	line(image, rects.at(i).at(0), rects.at(i).at(2), Scalar(0,255,0), 1, 8, 0);
    	line(image, rects.at(i).at(1), rects.at(i).at(3), Scalar(0,255,0), 1, 8, 0);
    	line(image, rects.at(i).at(2), rects.at(i).at(3), Scalar(0,255,0), 1, 8, 0);
    }*/
    //cout << rect.at(0) << ";" << rect.at(1) << ";" << rect.at(2) << ";" << rect.at(3) << endl;
    /*namedWindow("Projet2-Detection",CV_WINDOW_AUTOSIZE);//fenetre pour afficher la couleur moyenne
    imshow("Projet2-Detection", img_object);
        
    waitKey(0);*/
    imwrite(name.append("temp.jpg"), image);
    return rect;
}

int nbNotColor(cv::Mat image, int i, int j, cv::Vec3b col)
{
	int res = 0;
	int marge = 30;
	if (i - 1 > 0 && !goodColor(image.at<cv::Vec3b>(i-1,j), col, marge))
		res = res + 1;
	if (i + 1 < image.rows && !goodColor(image.at<cv::Vec3b>(i+1,j), col, marge))
		res = res + 1;
	if (j - 1 > 0 && !goodColor(image.at<cv::Vec3b>(i,j-1), col, marge))
		res = res + 1;
	if (j + 1 < image.rows && !goodColor(image.at<cv::Vec3b>(i,j+1), col, marge))
		res = res + 1;
	return res;
}

vector<cv::Point2d> Zernike::getCanette(cv::Mat image, int marge, string name)
{
	cv::Vec3b col(18,19,164);
	vector<cv::Point2d> vect = detectObjectWithColor(image, col, marge, name);
	return vect;
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

double zernike(cv::Mat image, int n, int m)
{
	 double rho, theta;
    cmplx zernike(0,0);
    double rows = image.rows;
    double cols = image.cols;

    for(int i=0;i<rows;i++)
    {
        //map y to unit circle
        for(int j=0;j<cols;j++)
        {
                   //map x, y  to unit ci rcle
            double y = -1 + ((2*(j+0.5))/cols);
            double x = -1 + ((2*(i+0.5))/rows);
            theta=0.0;
            rho = (double)sqrt(x*x+y*y);

            if(rho<=1.0)
            {
					if(x!=0)
						theta = atan(y/x);

					cv::Vec3b info = image.at<cv::Vec3b>(i,j);
					int val = (int) 0.1*info.val[0]+0.6*info.val[1]+0.3*info.val[2];
					cmplx pixel = val;
					//cout << "i = " << i << " , j = " << j << " , x = " << x << " , y = " << y << " , rho = " << rho << " , val = " << val << " , pixel = " << pixel << endl;
					zernike+= conj( R(n,m,rho)* cmplx(0, m*theta))*pixel; //Znm = [Rnm(r)*exp(-jmO]]*f(x,y)
            }
        }
    }
    cmplx val = (n+1)/M_PI;
    zernike= zernike*val;

    double result = abs(zernike/sqrt(image.rows*image.cols)); //return normalised Zernike moment   
    return result;
}
					
double Zernike::CalculateMoments(cv::Mat image)
{
	double res = 0;
    for(unsigned int n=0; n < 5; n++)
    	for(unsigned int m=0; m < 3;m++)
    	{
    		if ((n-m) % 2 == 0 && m <= n)
    		{
    			res = res + zernike(image, n, m);
    		}	
    	}
   return res;
}