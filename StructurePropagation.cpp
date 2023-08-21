#include "stdafx.h"
#include "StructurePropagation.h"
#include "OpenCvUtility.h"

int Dist(Vec3b x, Vec3b y)
{
	return (x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1]) + (x[2] - y[2]) * (x[2] - y[2]);
}

int Dist(pair<int, int> x, pair<int, int> y)
{
	return (x.first - y.first) * (x.first - y.first) + (x.second - y.second) * (x.second - y.second);
}

vector<Vec3b> GetEnvironment(int x, int y, Mat mat)
{
	vector<Vec3b> e1(4);
	x = MIN(mat.cols - 2, x);
	x = MAX(1, x);
	y = MIN(mat.rows - 2, y);
	y = MAX(1, y);

	e1[0] = mat.at<Vec3b>( y - 1,x-1);
	e1[1] = mat.at<Vec3b>(y - 1,x);
	e1[2] = mat.at<Vec3b>( y - 1,x+1);
	e1[3] = mat.at<Vec3b>( y,x-1);
	return e1;
}

vector<Vec3b> GetEnvironment2(int x, int y, Mat mat)
{
	vector<Vec3b> e1(4);
	e1[0] = mat.at<uchar>(y - 1, x - 1);
	e1[1] = mat.at<uchar>(y - 1, x);
	e1[2] = mat.at<uchar>(y - 1, x + 1);
	e1[3] = mat.at<uchar>(y, x - 1);
	return e1;
}

vector<Vec3b> GetEnvironment(pair<int, int> p[4], Mat mat)
{
	vector<Vec3b> v;
	for (int i = 0; i < 4; i++)
	{
		v.push_back(mat.at<Vec3b>(p[i].second, p[i].first));
	}
	return v;
}

int GetMinEnv(vector<Vec3b> x, vector<vector<Vec3b>> y)
{
	int minV = INT_MAX, minI = 0, temp;
	for (int i = 0; i < y.size(); i++)
	{
		temp = 0;
		if (y[i].size() == 0) continue;
		for (int j = 0; j < x.size(); j++)
		{
			temp += Dist(x[j], y[i][j]);
		}
		if (temp < minV)
		{
			minV = temp;
			minI = i;
		}
	}
	return minI;
}

void StructurePropagation::SetParm(int _blocksize, int _samplestep, int _iscurve)
{
	this->blocksize = _blocksize;
	this->samplestep = _samplestep;
	this->iscurve = _iscurve;
}

bool StructurePropagation::InsideMask(Point p)
{
	if (mask.at<uchar>(p) != 0) return false;
	else return true;
}

bool StructurePropagation::InsideMask(int x, int y)
{
	if (mask.at<uchar>(y,x) == 0) return true;
	else return false;
}

bool StructurePropagation::InsidePatch(Patch patch, Point p)
{
	if (p.x<patch.Location.x - blocksize / 2 || p.x>patch.Location.x + blocksize / 2) return false;
	if (p.y<patch.Location.y - blocksize / 2 || p.y>patch.Location.y + blocksize / 2) return false;
	return true;
}

bool StructurePropagation::OnMaskBorder(Point p)
{
	int px=p.x, py=p.y;
	int hasout=0, hasin=0;
	for (int x = MAX(px - blocksize / 2,0); x < MIN(px + blocksize / 2,mat.cols); x++)
	{
		for (int y = MAX(py - blocksize / 2,0); y < MIN(py + blocksize / 2,mat.rows); y++)
		{
			if (x == px && y == py) continue;
			if (InsideMask(x, y)) hasin = 1;
			else hasout = 1;
		}
	}
	if (hasin * hasout == 1) return true;
	else return false;
}

void StructurePropagation::Sampling(const vector<Point>& points)
{
	int index = 0,samplecnt=1,targetcnt=1;
	const int targetstep = blocksize / 2;
	for (int i=0;i<points.size();i++)
	{
		auto it = points[i];
		if (it.x > mat.cols || it.x < 0) continue;
		if (it.y > mat.rows || it.y < 0) continue;
		if (OnMaskBorder(it))
		{
			
			if (targetcnt == 1)
			{
				Targets.push_back(Patch(it, index++, Patch::ON_BORDER));
				targetcnt++;
				vector<Point> line;
				for (int j = 0; j < points.size(); j++)
				{
					if (InsidePatch(Targets.back(), points[j])) line.push_back(points[j]);
					else if (j > i + 1) break;
				}
				Targets.back().line = line;
			}
			else targetcnt = (targetcnt + 1) % targetstep;
		}
		else if (InsideMask(it))
		{
			if (targetcnt == 1)
			{
				Targets.push_back(Patch(it, index++, Patch::INSIDE_BORDER));
				targetcnt++;
				vector<Point> line;
				for (int j = 0; j < points.size(); j++)
				{
					if (InsidePatch(Targets.back(), points[j])) line.push_back(points[j]);
					else if (j > i + 1) break;
				}
				Targets.back().line = line;
			}
			else targetcnt = (targetcnt + 1) % targetstep;
		}
		else
		{
			if (samplecnt == 1)
			{
				Samples.push_back(Patch(it, index++, Patch::OUTSIDE_BORDER));
				samplecnt++;
				vector<Point> line;
				for (int j = 0; j < points.size(); j++)
				{
					if (InsidePatch(Samples.back(), points[j])) line.push_back(points[j]);
					else if (j > i + 1) break;
				}
				Samples.back().line = line;
			}
			else samplecnt = (samplecnt + 1) % samplestep;
		}
	}

}

void StructurePropagation::BsSampling(const vector<Point>& points)
{
	int lastx=-1, lasty=-1;
	int index = 0, samplecnt = 1, targetcnt = 1;
	Point lastp=points[0];
	BPoint lastbp;
	Point p;
	int isbeginpoint = 0;
	int targetstep = blocksize / 2;
	cout << "points size : " << points.size() << endl;
	cout << "targetstep is " << targetstep << endl;
	cout << "samplestep is " << samplestep << endl;
	for (int i = 0; i < points.size(); i++)
	{
		isbeginpoint = 0;
		p = points[i];
		if (i == points.size() - 1)
		{
			if (points[i].y < 0)
				EndPoints.push_back(i - 1);
			else
				EndPoints.push_back(i);
		}
		if (p.x > mat.cols || p.x < 0) continue;
		if (p.y > mat.rows || p.y < 0) continue;
		if (i > 0)
		{
			if (p.x == lastx && p.y == lasty) continue;
		}
		if (i == 0) isbeginpoint = 1;
		else if ((p.x - lastx) * (p.x - lastx) + (p.y - lasty) * (p.y - lasty) > 2) isbeginpoint = 1;
		if(i>0)cout << "( " << p.x << "," << p.y << ")   Dist = "<< (p.x - lastx) * (p.x - lastx) + (p.y - lasty) * (p.y - lasty) << endl;
		lastx = p.x;
		lasty = p.y;


		if (isbeginpoint) {
			targetcnt = 1; BeginPoints.push_back(i);
		}
		if ((isbeginpoint && i > 0))
		{
			EndPoints.push_back(i - 1);
		}


		if (isbeginpoint) cout << "new begin" << endl;

		if (OnMaskBorder(p))
		{
			if (targetcnt == 1)
			{
				Targets.push_back(Patch(p, index++, Patch::ON_BORDER));
				targetcnt++;
				vector<Point> line;
				for (int j = 0; j < points.size(); j++)
				{
					if (InsidePatch(Targets.back(), points[j])) line.push_back(points[j]);
					else if (j > i + 1) break;
				}
				Targets.back().line = line;

				BPoint bp(p);
				bp.TargetIndex = Targets.size();
				if (i > 0)
				{
					if(bp.AddNeighbor(lastbp,blocksize)) lastbp.AddNeighbor(bp,blocksize);
				}
				if (isbeginpoint || pointmap[pair<int, int>(p.x, p.y)] > 0)
				{
					for (int i = 0; i < Bpoints.size(); i++)
					{
						if (bp.AddNeighbor(Bpoints[i],blocksize)) Bpoints[i].AddNeighbor(bp,blocksize);
					}
				}
				Bpoints.push_back(bp);
				pointmap[pair<int, int>(p.x, p.y)]++;
				lastbp = bp;
				cout << "points[" << i << "] Neighbors: " << Bpoints.back().Neighbors.size() << endl;
			}
			else targetcnt = (targetcnt + 1) % targetstep;
		}
		else if (InsideMask(p))
		{
			if (targetcnt == 1)
			{
				Targets.push_back(Patch(p, index++, Patch::INSIDE_BORDER));
				targetcnt++;
				vector<Point> line;
				for (int j = 0; j < points.size(); j++)
				{
					if (InsidePatch(Targets.back(), points[j])) line.push_back(points[j]);
					else if (j > i + 1) break;
				}
				Targets.back().line = line;

				BPoint bp(p);
				bp.TargetIndex = Targets.size();
				if (i > 0)
				{
					if (bp.AddNeighbor(lastbp,blocksize)) lastbp.AddNeighbor(bp,blocksize);
				}
				if (isbeginpoint|| pointmap[pair<int, int>(p.x, p.y)] > 0)
				{
					for (int i = 0; i < Bpoints.size(); i++)
					{
						if (bp.AddNeighbor(Bpoints[i],blocksize)) Bpoints[i].AddNeighbor(bp,blocksize);
					}
				}
				Bpoints.push_back(bp);
				pointmap[pair<int, int>(p.x, p.y)]++;
				lastbp = bp;
				cout << "points[" << i << "] Neighbors: " << Bpoints.back().Neighbors.size() << endl;
			}
			else targetcnt = (targetcnt + 1) % targetstep;
		}
		else
		{
			if (samplecnt == 1)
			{
				Samples.push_back(Patch(p, index++, Patch::OUTSIDE_BORDER));
				samplecnt++;
				vector<Point> line;
				for (int j = 0; j < points.size(); j++)
				{
					if (InsidePatch(Samples.back(), points[j])) line.push_back(points[j]);
					else if (j > i + 1) break;
				}
				Samples.back().line = line;
			}
			else samplecnt = (samplecnt + 1) % samplestep;
		}
	}
	cout << "Bpoints size is " << Bpoints.size() << endl;
}

int StructurePropagation::PatchSSD(Point xi1,Point xi2,Point offset)
{
	int x, y, xp, yp;
	int ssd=0, cnt = 0;
	for (y = -blocksize / 2; y < blocksize / 2; y++)
	{
		yp = y + offset.y;
		if (yp<-blocksize / 2 || yp>blocksize / 2) continue;
		const auto* ptr1 = mat.ptr<Vec3b>(y + xi1.y);
		const auto* ptr2 = mat.ptr<Vec3b>(yp + xi2.y);
		for (x = -blocksize / 2; x < blocksize / 2; x++)
		{
			xp = x + offset.x;
			if (xp<-blocksize / 2 || xp>blocksize / 2) continue;
			ssd += Dist(ptr1[xi1.x + x], ptr2[xi2.x + xp]);
			cnt++;
		}
	}
	if (cnt == 0) return 0;
	return ssd / cnt;
}

int StructurePropagation::PatchSSD(Point i, Point xi)
{
	int ssd = 0;
	int cnt = 0;
	for (int y = -blocksize / 2; y < blocksize / 2; y++)
	{
		const auto* ptri = mat.ptr<Vec3b>(y + i.y);
		const auto* ptrxi = mat.ptr<Vec3b>(y + xi.y);
		for (int x = -blocksize / 2; x < blocksize / 2; x++)
		{
			if (ptri[x + i.x] != Vec3b(0,0,0))
			{
				ssd += Dist(ptri[x + i.x], ptrxi[x + xi.x]);
				cnt++;
			}
		}
	}
	if (cnt == 0) return 0;
	return ssd / cnt;
}

int StructurePropagation::Es(int S, int T)
{
	Patch spatch = Samples[S];
	Patch tpatch = Targets[T];

	Point spoint = spatch.Location;
	Point tpoint = tpatch.Location;

	vector<Point> sline = spatch.line;
	vector<Point> tline = tpatch.line;

	int result=0;
	int distance;
	int temp;

	for (int i = 0; i < sline.size(); i++)
	{
		distance = INT_MAX;
		for (int j = 0; j < tline.size(); j++)
		{
			temp = (sline[i].x - tline[j].x) * (sline[i].x - tline[j].x) + (sline[i].y - tline[j].y) * (sline[i].y - tline[j].y);
			if (temp < distance) distance = temp;
		}
		result += distance;
	}
	for (int i = 0; i < tline.size(); i++)
	{
		distance = INT_MAX;
		for (int j = 0; j < sline.size(); j++)
		{
			temp = (tline[i].x - sline[j].x) * (tline[i].x - sline[j].x) + (tline[i].y - sline[j].y) * (tline[i].y - sline[j].y);
			if (temp < distance) distance = temp;
		}
		result += distance;
	}
	return result / (tline.size() + sline.size());
}

int StructurePropagation::Ei(int S, int T)
{
	Patch spatch = Samples[S];
	Patch tpatch = Targets[T];
	
	Point spoint = spatch.Location;
	Point tpoint = tpatch.Location;

	if (tpatch.border == tpatch.ON_BORDER)
	{
		return PatchSSD(tpoint, spoint);
	}
	else return 0;
}

int StructurePropagation::E1(int S, int T)
{
	return int(0.25 * Es(S, T) + 0.75 * Ei(S, T));
}

int StructurePropagation::E2(int S1,int S2,int T1, int T2)
{
	Patch s1patch = Samples[S1];
	Patch s2patch = Samples[S2];

	Point s1point = s1patch.Location;
	Point s2point = s2patch.Location;

	Patch t1patch = Targets[T1];
	Patch t2patch = Targets[T2];

	Point t1point = t1patch.Location;
	Point t2point = t2patch.Location;
	Point offset = t2point - t1point;

	return PatchSSD(s1point, s2point, offset);
}

void StructurePropagation::DP()
{
	int L = Samples.size();
	int H = Targets.size();
	int* M = new int[L * H];
	int* record = new int[L * H];
	for (int i = 0; i < H; i++)
	{
		M[i] = E1(0, i);
	}
	int minD;
	int D = 0;
	cout << "DP levels is " << H << " : " << endl;
	for (int i = 1; i < H; i++)
	{
		cout << "DP: level " << i << " done..." << endl;
		for (int j = 0; j < L; j++)
		{
			minD = INT_MAX;
			for (int k = 0; k < L; k++)
			{
				D = E2(k, j, i - 1, i) + M[(i - 1) * L + k];
				if (D < minD) 
				{ 
					minD = D;
					record[i * L + j] = k;
				}
				else continue;
			}
			M[i * L + j] = E1(j, i) + minD;
		}
	}
	vector<int> DP_temp;
	int minindex = 0,minM=INT_MAX;
	for (int i = 0; i < L; i++)
	{
		if (M[(H - 1) * L + i] < minM)
		{
			minM = M[(H - 1) * L + i];
			minindex = i;
		}
	}
	for (int i = H-1; i >= 0; i--)
	{
		DP_temp.push_back(minindex);
		minindex = record[i * L + minindex];
	}
	for (int i = 0; i < H; i++)
	{
		DP_result.push_back(DP_temp.back());
		DP_temp.pop_back();
	}
}

void StructurePropagation::BSDP()
{
	vector<vector<int>> M,M_c;
	vector<pair<int, int>> pairs;
	map<pair<int, int>,int> toindex;
	vector<BPoint> Ns;
	int Sample_size = Samples.size();
	int index = 0;
	for (int i = 0; i < Bpoints.size(); i++)
	{
		Ns = Bpoints[i].GetNeighbors();
		for (int j = 0; j < Ns.size(); j++)
		{
			vector<int> temp;
			for (int k = 0; k < Sample_size; k++)
			{
				temp.push_back(0);
			}
			M.push_back(temp);
			M_c.push_back(temp);

			pairs.push_back(pair<int, int>(i, Ns[j].TargetIndex));
			toindex[pair<int, int>(i, Ns[j].TargetIndex)] = index++;
		}
	}
	cout << "initialize done" << endl;
	int fromp, top;
	int temp, minE1, minM;
	int* TargetRecord = new int[Bpoints.size()];
	int cnt;
	cout << "step 2 begin" << endl;
	for (int t = 1; t < Bpoints.size(); t++)
	{
		cout << "iteration " << t << " begin" << endl;
		cout << "DP for M levels: "<<pairs.size() << endl;
		for (int i = 0; i < pairs.size(); i++)
		{
			cout << "DP on " << i << " level" << endl;
			fromp = pairs[i].first;
			top = pairs[i].second;
			if (t == 1)
			{
				minE1 = INT_MAX;
				for (int j = 0; j < Samples.size(); j++)
				{
					temp = E1(j, fromp);
					if (temp < minE1)
					{
						minE1 = temp;
						TargetRecord[fromp] = j;
					}
				}
				for (int j = 0; j < Sample_size; j++)
				{
					M[i][j] = minE1;
				}
			}
			else
			{
				for (int j = 0; j < Sample_size; j++)
				{
					minM = INT_MAX;
					for (int k = 0; k < Sample_size; k++)
					{
						temp = E1(k, fromp) + E2(k, j, fromp, top);
						for (int kk = 0; kk < Bpoints[top].Neighbors.size(); kk++)
						{
							if (Bpoints[top].Neighbors[kk].TargetIndex != j)
							{
								temp += M_c[toindex[pair<int, int>(kk, i)]][k];
							}
						}
						if (temp < minM) minM = temp;
					}
					M[i][j] = minM;
				}
			}
		}
		cnt = 0;
		for (int i = 0; i < M.size(); i++)
		{
			for (int j = 0; j < M[i].size(); j++)
			{
				if (M[i][j] != M_c[i][j])
				{
					M_c[i][j] = M[i][j];
					cnt++;
				}
			}
		}
		cout << "cnt = " << cnt << endl;
		if (cnt == 0) break;
	}
	int minx = 0, minV = INT_MAX;
	for (int i = 0; i < Targets.size(); i++)
	{
		minx = 0;
		minV = INT_MAX;
		for (int j = 0; j < Sample_size; j++)
		{
			temp = E1(j, i);
			for (int k = 0; k < Bpoints[i].Neighbors.size(); k++)
			{
				temp += M[toindex[pair<int, int>(Bpoints[i].Neighbors[k].TargetIndex, i)]][j];
			}
			if (temp < minV)
			{
				minV = temp;
				minx = j;
			}
		}
		DP_result.push_back(minx);
	}
}

int StructurePropagation::GetMaxDist(int begin,vector<int>& mark)
{
	int Dist = 0;
	int LeftMaxDist = 0, temp;
	vector<int> record(mark);
	BPoint bp = Bpoints[begin];

	while (true)
	{
		mark[bp.TargetIndex] = 1;
		if (bp.Neighbors.size() > 2)
		{
			LeftMaxDist = 0;
			for (int j = 0; j < bp.Neighbors.size(); j++)
			{
				if (mark[bp.Neighbors[j].TargetIndex] == 0)
				{
					temp = GetMaxDist(bp.Neighbors[j].TargetIndex, mark);
					if (temp > LeftMaxDist) LeftMaxDist = temp;
				}
			}
			Dist += LeftMaxDist;
			break;
		}
		else if(bp.Neighbors.size()==1)
		{
			Dist++;
			break;
		}
		else
		{
			Dist++;
			if (mark[bp.Neighbors[0].TargetIndex] == 0) bp = bp.Neighbors[0];
			else bp = bp.Neighbors[1];
		}
	}
	for (int i = 0; i < mark.size(); i++)
	{
		mark[i] = record[i];
	}
	return Dist;
}

void StructurePropagation::BsTest()
{
	int index = 0;
	for (int i = 15; i < 245; i++)
	{
		Bpoints.push_back(BPoint(Point(15, i)));
		Bpoints.back().TargetIndex = index++;
	}
	for (int i = 15; i < 245; i++)
	{
		Bpoints.push_back(BPoint(Point(45, i)));
		Bpoints.back().TargetIndex = index++;
	}
	for (int i = 10; i < 60; i++)
	{
		Bpoints.push_back(BPoint(Point(i, 50)));
		Bpoints.back().TargetIndex = index++;
	}
	vector<int> mark;
	for (int i = 0; i < Bpoints.size(); i++)
	{
		mark.push_back(0);
	}
	cout << "Distance is " << GetMaxDist(0, mark) << endl;;
	
}

void StructurePropagation::CompletePatch(Mat& result, const Patch& tpatch,const Patch& spatch)
{
	const int tx = tpatch.Location.x;
	const int ty = tpatch.Location.y;

	const int sx = spatch.Location.x;
	const int sy = spatch.Location.y;

	int partial[4*3];
	Vec3b x1, x2, x0,value;
	int cnt;
	int L1, L2, D = blocksize / 2;

	for (int y = - blocksize / 2; y < blocksize / 2; y++)
	{
		auto* ptrt = result.ptr<Vec3b>(MIN(ty + y,result.rows-1));
		auto* ptrs = result.ptr<Vec3b>(MIN(sy + y,result.rows-1));
		for (int x = - blocksize / 2; x < blocksize / 2; x++)
		{
			if (ptrt[x + tx] != Vec3b(0, 0, 0)) ptrt[x + tx] = ptrt[x + tx] / 2 + ptrs[x + sx] / 2;
			else ptrt[x + tx] = ptrs[x + sx];
			//ptrt[x + tx] = ptrs[x + sx];
			//src[pair<int, int>(x + tx, y + ty)] = pair<int, int>(x + sx, y + sy);
		}
	}
}

void StructurePropagation::Completing(Mat& result)
{
	for (int i = 0; i < Targets.size(); i++)
	{
		CompletePatch(result, Targets[i], Samples[DP_result[i]]);
		cout << "Target " << i << " completed..." << endl;
	}
}

void StructurePropagation::TextureComplete(Mat& result, const Mat& _mask,const vector<Point> points)
{
	//Mat test1(result.size(), result.type());

	int number = 1;
	int begin, end, cnt;
	int dx, dy;
	int x, y;
	int lastx=-1, lasty=-1;
	const int rows = result.rows;
	const int cols = result.cols;
	Mat Number(result.size(),CV_8UC1,Scalar(0));
	for (int x = 0; x < cols; x++)
	{
		for (int y = 0; y < rows; y++)
		{
			if (InsideMask(x, y)) {
				Number.at<uchar>(y, x) = 1; // std::cout << "(" << x << "," << y << "); ";
			}
			else if (src.count(pair<int, int>(x, y)) == 0) src[pair<int, int>(x, y)] = pair<int, int>(x, y);
		}
	}
	for (int i = 0; i < BeginPoints.size(); i++)
	{
		begin = BeginPoints[i];
		end = EndPoints[i];
		dx = points[end].x - points[begin].x;
		dy = points[end].y - points[begin].y;
		cnt = 0;
		if (abs(dx) > abs(dy))
		{
			for (int j = begin; j <= end; j++)
			{
				if (points[j].y < 0) continue;
				x = points[j].x;
				if (j > begin && dx * (x - lastx) <= 0) continue;
				cnt = 0;
				for (y = points[j].y+blocksize/2;; y++)
				{
					if (y > rows - 2) break;
					Number.at<uchar>(y, x) += number;
					//cout << "(" << x << "," << y << "); ";
					//if (!InsideMask(x, y))cnt++;
				}
				for (y = points[j].y - blocksize / 2;; y--)
				{
					if (y <3) break;
					Number.at<uchar>(y, x) += number*2;
				}
				lastx = x;
			}

		}
		else
		{
			for (int j = begin; j <= end; j++)
			{
				if (points[j].y < 0) continue;
				y = points[j].y;
				if (j > begin && dy * (y - lasty) <= 0) continue;
				cnt = 0;
				for (x = points[j].x+blocksize/2;; x++)
				{
					if (cnt > 10 || x > cols - 2) break;
					Number.at<uchar>(y, x) += number;
					//if (!InsideMask(x, y))cnt++;
				}
				for (x = points[j].x - blocksize / 2;; x--)
				{
					if (x<3) break;
					Number.at<uchar>(y, x) += number*2;
				}
				lasty = y;
			}
		}
		number *= 4;

	}


	//Mat test(Number.size(), Number.type());
	//for (int y = 0; y < test.rows; y++)
	//{
	//	for (int x = 0; x < test.cols; x++)
	//	{
	//		test.at<uchar>(y, x) = Number.at<uchar>(y, x) * 16;
	//		int i1 = Number.at<uchar>(y, x);
	//		if (y == 400 && x == 100) cout <<"(100,400) = "<< i1 << endl;
	//		if (y == 100 && x == 100) cout << "(100,100) = " << i1 << endl;
	//	}
	//}
	//imshow("test", test);
	//waitKey(0);
	//cout << Number << endl;



	vector<vector<pair<int,int>>> Num2TP(number+1);
	vector<vector<pair<int,int>>> Num2SP(number+1);
	int* SPBound = new int[(number + 1)*4];
	for (int i = 0; i < number + 1; i++)
	{
		SPBound[i * 4 + 0] = INT_MAX;
		SPBound[i * 4 + 1] = 0;
		SPBound[i * 4 + 2] = INT_MAX;
		SPBound[i * 4 + 3] = 0;
	}
	int temp;
	int index = 0;
	for (int y = 0; y < rows; y++)
	{
		for (int x = 0; x < cols; x++)
		{
			temp = Number.at<uchar>(y, x);
			if (mask.at<uchar>(y, x) > 0)
			{
				Num2SP[temp].push_back(pair<int, int>(x, y));		
				if (x < SPBound[temp * 4 + 0]) SPBound[temp * 4 + 0] = x;
				if (x > SPBound[temp * 4 + 1]) SPBound[temp * 4 + 1] = x;
				if (y < SPBound[temp * 4 + 2]) SPBound[temp * 4 + 2] = y;
				if (y > SPBound[temp * 4 + 3]) SPBound[temp * 4 + 3] = y;
			}
			else
			{
				Num2TP[temp].push_back(pair<int, int>(x, y));
			}
		}
	}
	for (int i = 0; i < number + 1; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cout << SPBound[i * 4 + j] << " ";
		}
		cout << endl;
	}
	/*
	int rand[10] = { 1,3,9,8,14,5,72,66,22,29 };
	for (int n = 1; n <= number; n++)
	{
		int size = Num2SP[n-1].size();
		for (int i = 0; i < Num2TP[n].size(); i++)
		{
			src[Num2TP[n][i]] = Num2SP[n-1][(rand[(i * 23 + 9) % 10] * 17 + 37) % size];
			//src.insert(pair<pair<int,int>,pair<int,int>>(Num2TP[n][i],Num2SP[n][(rand[(i * 23 + 9) % 10] * 17 + 37)%size]));
		}
		for (int i = 0; i < Num2TP[n].size(); i++)
		{
			vector<Vec3b> e1(4);
			vector<vector<Vec3b>> e2(4);
			int x, y;
			x = Num2TP[n][i].first;
			y = Num2TP[n][i].second;
			pair<int, int> pairs[4];
			pairs[0] = pair<int, int>(x - 1, y - 1);
			pairs[1] = pair<int, int>(x, y - 1);
			pairs[2] = pair<int, int>(x + 1, y - 1);
			pairs[3] = pair<int, int>(x - 1, y);

			e1 = GetEnvironment(x, y, mat);
			for (int k = 0; k < 4; k++)
			{
				e2[k] = GetEnvironment(src[pairs[k]].first, src[pairs[k]].second, mat);
			}
			int index = GetMinEnv(e1, e2);
			int newx, newy;
			switch (index) {
			case 0:
				newx = src[pairs[0]].first + 1;
				newy = src[pairs[0]].second + 1;
				break;
			case 1:
				newx = src[pairs[1]].first;
				newy = src[pairs[1]].second + 1;
				break;
			case 2:
				newx = src[pairs[2]].first - 1;
				newy = src[pairs[2]].second + 1;
				break;
			case 3:
				newx = src[pairs[3]].first + 1;
				newy = src[pairs[3]].second;
				break;
			}
			newx = (newx >= 0) ? newx : 0;
			newx = (newx < cols) ? newx : cols - 1;
			newy = (newy >= 0) ? newy : 0;
			newy = (newy < rows) ? newy : rows - 1;
			src[pair<int, int>(x, y)] = pair<int, int>(newx, newy);
		    result.at<Vec3b>(y, x) = result.at<Vec3b>(newy, newx);
		}
	}
	*/

	/*
	int _maskrows = _mask.rows;
	int _maskcols = _mask.cols;
	cout << "initialize mask levels: " << _maskrows << endl;
	for (int y = 0; y < _maskrows; y++)
	{
		cout << "level : " << y << endl;
		for (int x = 0; x < _maskcols; x++)
		{
			if (_mask.at<uchar>(y, x) == 255)
			{
				if (_mask.at<uchar>(y, MIN(x + 1,_maskcols-1)) == 0 || _mask.at<uchar>(MIN(y + 1, _maskrows-1), MAX(x - 1,0)) == 0 ||
					_mask.at<uchar>(MIN(y + 1,_maskrows-1), x) == 0 || _mask.at<uchar>(MIN(y + 1, _maskrows-1), MIN(x + 1, _maskcols-1)) == 0)
				{
					int minxx = 0, minyy = 0, minV = INT_MAX, temp,breakflag=0;
					Vec3b v = mat.at<Vec3b>(y, x);
					for (int yy = y; yy >= 0 ; yy--)
					{
						breakflag = 0;
						for (int xx = (yy==y)?x-1:_maskcols; xx >= 0; xx--)
						{
							if (src.count(pair<int, int>(xx, yy))==0||src[pair<int,int>(xx,yy)]!=pair<int,int>(xx,yy)) continue;
							temp = Dist(mat.at<Vec3b>(yy, xx), v);
							if (temp == 0)
							{
								cout << "exactly equal!" << endl;
								minV = temp;
								minxx = xx;
								minyy = yy;
								breakflag = 1;
								break;
							}
							if (temp < minV)
							{
								minV = temp;
								minxx = xx;
								minyy = yy;
							}
						}
						if (breakflag == 1) break;
					}
					src[pair<int, int>(x, y)] = pair<int, int>(minxx, minyy);
				}
				else src[pair<int, int>(x, y)] = pair<int, int>(x, y);

			}
		}
	}

	for (int n = 1; n <= number; n++)
	{
		int size = Num2SP[n - 1].size();
		pair<int, int> x0, x1, x2;
		int rand1, rand2;
		for (int i = 0; i < Num2TP[n].size(); i++)
		{
			x0 = Num2TP[n][i];
			rand1 = rand() % size;
			rand2 = rand() % size;
			x1 = Num2SP[n - 1][rand1];
			x2 = Num2SP[n - 1][rand2];
			if (Dist(x1, x0) < Dist(x2, x0)) src[Num2TP[n][i]] = Num2SP[n - 1][rand1];
			else src[Num2TP[n][i]] = Num2SP[n - 1][rand2];
		}

		for (int i = 0; i < Num2TP[n].size(); i++)
		{
			vector<Vec3b> e1(4);
			vector<vector<Vec3b>> e2(4);
			int x, y;
			x = Num2TP[n][i].first;
			y = Num2TP[n][i].second;
			pair<int, int> pairs[4];
			x = MAX(x, 1);
			x = MIN(x, _maskcols - 2);
			y = MAX(y, 1);
			y = MIN(y, _maskrows - 2);
			pairs[0] = pair<int, int>(x - 1, y - 1);
			pairs[1] = pair<int, int>(x, y - 1);
			pairs[2] = pair<int, int>(x + 1, y - 1);
			pairs[3] = pair<int, int>(x - 1, y);
			e1 = GetEnvironment(x, y, mat);
			vector<int> temp_index;
			if (_mask.at<uchar>( y - 1,x - 1) > 0 && src.count(pairs[0]) == 1 && src[pairs[0]] != pairs[0])
				temp_index.push_back(0);
			if (_mask.at<uchar>( y - 1,x) > 0 && src.count(pairs[1]) == 1 && src[pairs[1]] != pairs[1])
				temp_index.push_back(1);
			if (_mask.at<uchar>( y - 1,x + 1) > 0 && src.count(pairs[2]) == 1 && src[pairs[2]] != pairs[2])
				temp_index.push_back(2);
			if (_mask.at<uchar>( y,x - 1) > 0 && src.count(pairs[3]) == 1 && src[pairs[3]] != pairs[3])
				temp_index.push_back(3);

			if (temp_index.size() > 0)
			{
				for (int k = 0; k < temp_index.size(); k++)
				{
					e2[temp_index[k]]= GetEnvironment(src[pairs[temp_index[k]]].first, src[pairs[temp_index[k]]].second, mat);
				}
			}
			else
			{
				for (int k = 0; k < 4; k++)
				{
					e2[k] = GetEnvironment(src[pairs[k]].first, src[pairs[k]].second, mat);
				}
			}
			int index = GetMinEnv(e1, e2);
			int newx, newy;
			switch (index) {
			case 0:
				newx = src[pairs[0]].first + 1;
				newy = src[pairs[0]].second + 1;
				break;
			case 1:
				newx = src[pairs[1]].first;
				newy = src[pairs[1]].second + 1;
				break;
			case 2:
				newx = src[pairs[2]].first - 1;
				newy = src[pairs[2]].second + 1;
				break;
			case 3:
				newx = src[pairs[3]].first + 1;
				newy = src[pairs[3]].second;
				break;
			}
			newx = (newx >= 0) ? newx : 0;
			newx = (newx < cols) ? newx : cols - 1;
			newy = (newy >= 0) ? newy : 0;
			newy = (newy < rows) ? newy : rows - 1;
			src[pair<int, int>(x, y)] = pair<int, int>(newx, newy);
			result.at<Vec3b>(y, x) = result.at<Vec3b>(newy, newx);
		}
	}
	*/
    int _maskcols = _mask.cols;
	int _maskrows = _mask.rows;
	for (int n = 1; n <= number; n++)
	{
		if (n == 7)
		{
			int a = 10;
			a++;
		}
		int size = Num2SP[n - 1].size();
		pair<int, int> x0, x1, x2;
		int rand1, rand2;
		for (int i = 0; i < Num2TP[n].size(); i++)
		{
			//x0 = Num2TP[n][i];
			//rand1 = rand() % size;
			//rand2 = rand() % size;
			//x1 = Num2SP[n - 1][rand1];
			//x2 = Num2SP[n - 1][rand2];
			//if (Dist(x1, x0) < Dist(x2, x0)) src[Num2TP[n][i]] = Num2SP[n - 1][rand1];
			//else src[Num2TP[n][i]] = Num2SP[n - 1][rand2];
			//src[Num2TP[n][i]] = Num2SP[n - 1][rand()%size];
		}

		for (int i = 0; i < Num2TP[n].size(); i++)
		{
			vector<Vec3b> e1(4);
			vector<vector<Vec3b>> e2(4);
			int x, y;
			x = Num2TP[n][i].first;
			y = Num2TP[n][i].second;
			pair<int, int> pairs[4];
			pair<int, int> srcs[4];
			x = MAX(x, 1);
			x = MIN(x, _maskcols - 2);
			y = MAX(y, 1);
			y = MIN(y, _maskrows - 2);
			pairs[0] = pair<int, int>(x - 1, y - 1);
			pairs[1] = pair<int, int>(x, y - 1);
			pairs[2] = pair<int, int>(x + 1, y - 1);
			pairs[3] = pair<int, int>(x - 1, y);
			float d,rand1;
			int tempx, tempy;
			for (int j = 0; j < 4; j++)
			{
				if (src.count(pairs[j]) == 0 || src[pairs[j]] == pairs[j])
				{
					srcs[j] = Num2SP[n - 1][rand() % size];
					tempx = srcs[j].first; tempy = srcs[j].second;
					d = ((float)(tempy - SPBound[(n - 1) * 4 + 2])) / (SPBound[(n - 1) * 4 + 3] - SPBound[(n - 1) * 4 + 2]);
					rand1 = pow(rand() % 100 / ((float)100), 4);
					while (rand1>d)
					{
						srcs[j] = Num2SP[n - 1][rand() % size];
						tempx = srcs[j].first; tempy = srcs[j].second;
						d = ((float)(tempy - SPBound[(n - 1) * 4 + 2])) / (SPBound[(n - 1) * 4 + 3] - SPBound[(n - 1) * 4 + 2]);
						rand1 = pow(rand() % 100 / ((float)100), 4);
					}
				}
				else srcs[j] = src[pairs[j]];
			}
			e1 = GetEnvironment(pairs,mat);

			tempx=srcs[0].first, tempy=srcs[0].second;
			while (tempx > cols - 3 || tempy > rows - 2||Number.at<uchar>(tempy+1,tempx+1)!=n-1)
			{
				srcs[0] = Num2SP[n - 1][rand() % size];
				tempx = srcs[0].first; tempy = srcs[0].second;
			}
			e2[0].push_back(mat.at<Vec3b>(tempy, tempx));
			e2[0].push_back(mat.at<Vec3b>(tempy, tempx + 1));
			e2[0].push_back(mat.at<Vec3b>(tempy, tempx + 2));
			e2[0].push_back(mat.at<Vec3b>(tempy + 1, tempx));

			tempx = srcs[1].first; tempy = srcs[1].second;
			while (tempx > cols - 2 || tempy > rows - 2 || tempx <1 || Number.at<uchar>(tempy + 1, tempx) != n - 1)
			{
				srcs[1] = Num2SP[n - 1][rand() % size];
				tempx = srcs[1].first; tempy = srcs[1].second;
			}
			e2[1].push_back(mat.at<Vec3b>(tempy, tempx - 1));
			e2[1].push_back(mat.at<Vec3b>(tempy, tempx));
			e2[1].push_back(mat.at<Vec3b>(tempy, tempx + 1));
			e2[1].push_back(mat.at<Vec3b>(tempy + 1, tempx - 1));

			tempx = srcs[2].first; tempy = srcs[2].second;
			while (tempy > rows - 2 || tempx < 2 || Number.at<uchar>(tempy + 1, tempx - 1) != n - 1)
			{
				srcs[2] = Num2SP[n - 1][rand() % size];
				tempx = srcs[2].first; tempy = srcs[2].second;
			}
			e2[2].push_back(mat.at<Vec3b>(tempy, tempx - 2));
			e2[2].push_back(mat.at<Vec3b>(tempy, tempx - 1));
			e2[2].push_back(mat.at<Vec3b>(tempy, tempx));
			e2[2].push_back(mat.at<Vec3b>(tempy + 1, tempx - 2));

			tempx = srcs[3].first; tempy = srcs[3].second;
			while (tempx > cols - 3 || tempy < 1 || Number.at<uchar>(tempy, tempx + 1) != n - 1)
			{
				srcs[3] = Num2SP[n - 1][rand() % size];
				tempx = srcs[3].first; tempy = srcs[3].second;
			}
			e2[3].push_back(mat.at<Vec3b>(tempy - 1, tempx));
			e2[3].push_back(mat.at<Vec3b>(tempy - 1, tempx + 1));
			e2[3].push_back(mat.at<Vec3b>(tempy - 1, tempx + 2));
			e2[3].push_back(mat.at<Vec3b>(tempy, tempx));


			//for (int k = 0; k < 4; k++)
			//{
			//	e2[k] = GetEnvironment(srcs[k].first, srcs[k].second, mat);
			//	//e2[k] = GetEnvironment(srcs[k], mat);
			//}
			int index = GetMinEnv(e1, e2);
			int newx, newy;
			switch (index) {
			case 0:
				newx = srcs[0].first + 1;
				newy = srcs[0].second + 1;
				break;
			case 1:
				newx = srcs[1].first;
				newy = srcs[1].second + 1;
				break;
			case 2:
				newx = srcs[2].first - 1;
				newy = srcs[2].second + 1;
				break;
			case 3:
				newx = srcs[3].first + 1;
				newy = srcs[3].second;
				break;
			}

			//if (Number.at<uchar>(newy, newx) != n-1)
			//{
			//	i--; continue;
			//}

			while (src[pair<int, int>(newx, newy)] != pair<int, int>(newx, newy))
			{
				newx = src[pair<int, int>(newx, newy)].first;
				newy = src[pair<int, int>(newx, newy)].second;
			}
			newx = (newx >= 0) ? newx : 0;
			newx = (newx < cols) ? newx : cols - 1;
			newy = (newy >= 0) ? newy : 0;
			newy = (newy < rows) ? newy : rows - 1;
			src[pair<int, int>(x, y)] = pair<int, int>(newx, newy);
			result.at<Vec3b>(y, x) = result.at<Vec3b>(newy, newx);
			//cout << "(" << newx << "," << newy << ") => (" << x << "," << y << ")" << endl;
		}
	}
}

void StructurePropagation::Run(const Mat1b& _mask, const Mat& _img, const vector<Point>& points, Mat& result)
{
	time_t time1, time2, time3,time4;
	time(&time1);
	this->mask = _mask;
	this->mat = _img;
	//BsSampling(points);
	Sampling(points);
	cout << "Sampling done..." << endl;
	DP();
	time(&time2);
	//BSDP();
	cout << "DP done..." << endl;
	Completing(result);
	cout << "structure compete done..." << endl;
	cout << "mask:" << endl;
	//cout << _mask << endl;
	for (int y = 0; y < mask.rows; y++)
	{
		for (int x = 0; x < mask.cols; x++)
		{
			if (_mask.at<uchar>(y, x) == 0 && result.at<Vec3b>(y, x) != Vec3b(0, 0, 0))
			{
				//cout << "(" << x << "," << y << "); ";
				mask.at<uchar>(y, x) = 255;
			}
		}
	}
	//imshow("runmask", mask);
	//cout << "\n";
	time(&time3);
	//TextureComplete(result, _mask, points);
	time(&time4);
	cout << "DP Time: " << difftime(time3, time2) << ";    Texture Time: " << difftime(time4, time3) << endl;
}