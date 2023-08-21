#pragma  once
#include "OpenCvUtility.h"

class Patch
{
public:
	enum BorderType { INVALID_BORDER = -1, OUTSIDE_BORDER = 0, ON_BORDER = 1, INSIDE_BORDER = 2 } ;
	Patch(Point _Location, int _index) :Location(_Location), index(_index) { border = INVALID_BORDER; }
	Patch(Point _Location, int _index, BorderType _type) :Location(_Location), index(_index), border(_type) { ; }
	Point Location;
	int index;
	BorderType border;
	vector<Point> line;
};

class BPoint : public Point 
{
public:
	BPoint(Point p) :Point(p) { TargetIndex = -1; }
	BPoint() :Point() { TargetIndex = -1; }
	int TargetIndex;
	vector<BPoint> Neighbors;
	vector<BPoint> GetNeighbors() { return Neighbors; }
	bool AddNeighbor(BPoint p,int maxdist)
	{
		if (p.x == this->x && p.y == this->y) return false;
		for (int i = 0; i < Neighbors.size(); i++)
		{
			if (p.x == Neighbors[i].x && p.y == Neighbors[i].y) return false;
		}
		if ((p.x - x) * (p.x - x) + (p.y - y) * (p.y - y) <= maxdist*maxdist)
		{
			Neighbors.push_back(p); cout << "add one neighbor" << endl; return true;
		}
		return false;
	}
};

class StructurePropagation
{
public:
	~StructurePropagation() {}
	void Run(const Mat1b& _mask, const Mat& _img, const vector<Point>& points, Mat& result);
	//利用目前类中已经存储的数据继续经行修补
	//void RunAgain();
	void SetParm(int _blocksize, int _samplestep, int _iscurve);
	void BsTest();
private:
	bool iscurve;
	const vector<Point>* Curve_Point;

	vector<BPoint> Bpoints;
	vector<Patch> Targets;
	vector<Patch> Samples;
	vector<int> BeginPoints;
	vector<int> EndPoints;
	map<pair<int, int>, int> pointmap;
	map<pair<int, int>, pair<int, int>> src;

	vector<int> DP_result;
	Mat1i DP_Mat;//col 是Target中的位置，row是Sample

	Mat1b mask;
	Mat mat;

	int blocksize;
	int samplestep;
	int GetMaxDist(int begin, vector<int>& mark);
	void Sampling(const vector<Point>& points);//目标区域采样以及修补块采样
	void BsSampling(const vector<Point>& points);


	int PatchSSD(Point i, Point xi);//块A与块B数据块以及相对位置
	int PatchSSD(Point i1,Point xi1,Point offset);

	void PointMapInit()
	{
		for (int x = 0; x < mat.cols; x++)
		{
			for (int y = 0; y < mat.rows; y++)
			{
				pointmap[pair<int, int>(x, y)] = 0;
			}
		}
	}

	int Es(int S, int T);
	int Ei(int S, int T);
	int E1(int S, int T);
	int E2(int S1,int S2,int T1, int T2);

	int Estructure(int S, int T);//T,S分别为Target与Sample的位置索引
	int Esmooth(Point pairA, Point pairB);//DP_Mat中的相邻列点，x为行，y为列

	void DP();
	void BSDP();
	void Completing(Mat& result);
	void CompletePatch(Mat& result, const Patch& tpatch, const Patch& spatch);
	void TextureComplete(Mat& result, const Mat& mask, const vector<Point> points);
	void GetBorder();

	bool InsideMask(Point p);
	bool InsideMask(int x, int y);
	bool OnMaskBorder(Point p);

	bool InsidePatch(Patch patch, Point p);

};