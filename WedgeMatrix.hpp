#ifndef WEDGEMATRIX_HPP
#define WEDGEMATRIX_HPP

#include <vector>

using namespace std;

class WedgeMatrix
{
public:
	WedgeMatrix(int k, bool full);
	virtual ~WedgeMatrix();
	
public:
	void ResolveFull();
	void ResolveParetos();
	void NormaliseStructs();
	void ReWeightStructs();
	void DeleteMemory();
	void SetNull(WedgeMatrix * nullMat);
	void SetAsNull();
	void MergeIn();
	bool IsNull();
	void Lock();
	void Unlock();
	bool IsLocked();
	
	void SubtractNull(bool doSD);
	
	void SendStructZ(short depth, vector<bool> * npath, float s, unsigned int freq);
	void SendParetoZ(short depth, vector<bool> * npath, float s, unsigned int freq);
	
	void SendPureStruct(short depth, float s, unsigned int freq);
	void SendPurePareto(short depth, float s, unsigned int freq);
	
	void TSendStructZ(short depth, short xn, short yn, short nn, float s, unsigned int freq);
	void TSendParetoZ(short depth, short xn, short yn, short nn, float s, unsigned int freq);
	void TSendBoth(short depth, short xn, short yn, short nn, float s, unsigned int freq);
	
	float GetScoreM(short i, short j, short k);
	float GetScoreD(short i, short j, short k);
	float GetZScoreM(short i);
	float GetZScoreD(short i);
	
	float GetNScoreS(short d, short n);
	float GetNScoreF(short d, short n);
	float GetNScoreShape(short d, short n);
	float GetNScorePN(short d, short n);
	
	unsigned int GetFreqM(short i, short j, short k);
	unsigned int GetDepthFreqM(short i, short j, short k);
	unsigned int GetShapeN(short i, short j, short k);
	
	double GetSumStruct();
	double GetNormalisedS();
	double GetSumP();
	double GetSumFs();
	double GetParetoN();
	
	void AddNFreq(unsigned int f);
	void AddHeadFreq(unsigned int f);
	
protected:
	int _k;
	bool _full;
	bool _hasNull;
	bool _isNull;
	bool _isLocked;
	
	double sumF;
	double sumS;
	double sumN;
	double sumP;
	
	unsigned int _NFReqs;
	unsigned int _headFreq;

	vector<float> **** _sMemory;
	vector<int> **** _fMemory;
	
	vector<float> ** _zFMemory;
	vector<int> ** _zNMemory;
	
	vector<float> _zStruct;
	vector<float> _zDev;
	
	vector<float> _zFTracker;
	vector<float> _zNTracker;
	
	float *** _sMeans;
	float *** _sDeviations;
	
	vector<float> ** _nMemoryS;
	vector<float> ** _nMemoryF;
	
	vector<float> ** _nMemoryShape;
	vector<float> ** _nMemoryPN;
	
	unsigned int *** _fTracker;
	unsigned int *** _dfTracker;
	unsigned int *** _nTracker;
	
	void FindMSD(vector<float> * smem, vector<int> * fmem, float * meanVal, float * sdVal);
	void MinZero(float * val);
	
	WedgeMatrix * _referenceNull;
	
	double CorrectFreq(unsigned int frq, short depth);

public:
	static float PLOIDY;
	static unsigned int SETSIZE;
	
	static void FindIndices(vector<bool> * npath, short * xn, short * yn, short * nn);
};

#endif // WEDGEMATRIX_HPP
