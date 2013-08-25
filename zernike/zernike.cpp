#include "zernike.h"
 
using namespace std; 
using namespace cv;
 
Zernike::Zernike()
{
} 

Zernike::~Zernike()
{
}

std::vector<cv::Point2d> getRect(std::vector<cv::Point2d> vect)
{
	vector<Point2d> res;
	double minX = -1;
	double minY = -1;
	double maxX = -1;
	double maxY = -1;
	for (unsigned int i=1; i < vect.size(); i++)
	{
		Point2d pt = vect.at(i);
		if (minX == -1 || minX > pt.x)
			minX = pt.x;
		if (minY == -1 || minY > pt.y)
			minY = pt.y;
		if (maxX == -1 || maxX < pt.x)
			maxX = pt.x;
		if (maxY == -1 || maxY < pt.y)
		maxY = pt.y;
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

void extractRect(Mat image, std::vector<cv::Point2d> rect)
{
	double minX = rect.at(0).x;
	double minY = rect.at(0).y;
	double maxX = rect.at(3).x;
	double maxY = rect.at(3).y; 
	Mat M(maxY-minY+1,maxX-minX+1, CV_8UC3);
	for (unsigned int i = minX; i <= maxX;i++)
		for (unsigned int j = minY; j <= maxY;j++)
			M.at<Vec3b>(j-minY,i-minX) = image.at<Vec3b>(j,i);
	 if (M.rows > M.cols)
    {
    	Mat tmp(maxX-minX+1,maxY-minY+1, CV_8UC3);
    	rotate(M,-90, tmp);
    	M = tmp;
    }
	imwrite("extracted.jpg", M);
}

vector<Point2d> detectObjectWithColor(Mat image, Vec3b scalColor)
{
    int margeR = 30;
    int margeG = 30;
    int margeB = 30;
 
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
        if(scalPix.val[0] > (scalColor.val[0] - margeR) && (scalPix.val[0] < (scalColor.val[0] + margeR)))
        {
            if(scalPix.val[1] > (scalColor.val[1] - margeG) && (scalPix.val[1] < (scalColor.val[1] + margeG)))
            {
                if(scalPix.val[2] > (scalColor.val[2] - margeB) && (scalPix.val[2] < (scalColor.val[2] + margeB)))
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
                }
            else
            {
                scalPainter.val[0] = 0;
                scalPainter.val[1] = 0;
                scalPainter.val[2] = 0;
            }
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
    vector<Point2d> rect = getRect(points);
    extractRect(image, rect);
    line(image, rect.at(0), rect.at(1), Scalar(0,255,0), 1, 8, 0);
    line(image, rect.at(0), rect.at(2), Scalar(0,255,0), 1, 8, 0);
    line(image, rect.at(1), rect.at(3), Scalar(0,255,0), 1, 8, 0);
    line(image, rect.at(2), rect.at(3), Scalar(0,255,0), 1, 8, 0);
    cout << rect.at(0) << ";" << rect.at(1) << ";" << rect.at(2) << ";" << rect.at(3) << endl;
    //namedWindow("Projet2-Detection",CV_WINDOW_AUTOSIZE);//fenetre pour afficher la couleur moyenne
    //imshow("Projet2-Detection", image);
    imwrite("temp.jpg", image);
    
    waitKey(0);
    return rect;
}

bool goodColor(cv::Vec3b pix, cv::Vec3b col)
{
	double bleu = pix.val[0];
	double vert = pix.val[1];
	double rouge = pix.val[2];
	if (bleu >= (col.val[0]-25) && bleu <= (col.val[0]+25) && vert >= (col.val[1]-25) && vert <= (col.val[1]+25) && rouge >= (col.val[2]-25) && rouge <= (col.val[2]+25))
		return true;
	return false;
}

int nbNotColor(cv::Mat image, int i, int j, cv::Vec3b col)
{
	int res = 0;
	if (i - 1 > 0 && !goodColor(image.at<cv::Vec3b>(i-1,j), col))
		res = res + 1;
	if (i + 1 < image.rows && !goodColor(image.at<cv::Vec3b>(i+1,j), col))
		res = res + 1;
	if (j - 1 > 0 && !goodColor(image.at<cv::Vec3b>(i,j-1), col))
		res = res + 1;
	if (j + 1 < image.rows && !goodColor(image.at<cv::Vec3b>(i,j+1), col))
		res = res + 1;
	return res;
}

vector<cv::Point2d> Zernike::getCanette(cv::Mat image)
{
	int orient = -1;
	cv::Vec3b col(18,19,165);
	vector<cv::Point2d> vect = detectObjectWithColor(image, col);
	/*for (int i=0; i < image.rows; i++)
	{
		for (int j=0; j < image.cols; j++)
		{
			cv::Vec3b info = image.at<cv::Vec3b>(i,j);
			if (goodColor(info, col))
			{
				 if (orient == -1)
				 {
				 	orient = nbNotColor(image, i, j, col);
				 	cout << orient << endl;
				 	vect.push_back(cv::Point2d(i,j));
				 }
				 else
				 {
				 	if (orient == nbNotColor(image, i, j, col))
				 		vect.push_back(cv::Point2d(i,j));
				 }
			}
		}
	}*/
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
					
double Zernike::CalculateMoments(cv::Mat image, int n, int m)
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
					zernike+= conj( R(n,m,rho)* polar(1.0, m*theta))*pixel; //Znm = [Rnm(r)*exp(-jmO]]*f(x,y)
            }
        }
    }
    cmplx val = (n+1)/M_PI;
    zernike= zernike*val;

    double result = abs(zernike/sqrt(image.rows*image.cols)); //return normalised Zernike moment
    return result;
}