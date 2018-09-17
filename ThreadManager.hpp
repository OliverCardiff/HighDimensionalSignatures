#ifndef THREADMANAGER_HPP
#define THREADMANAGER_HPP
#include <thread>
#include "WedgeMatrix.hpp"

class ThreadManager
{
public:
	
	struct Submission
	{
		unsigned int _freq;
		short _depth;
		double _s;
		short _xn;
		short _yn;
		short _nn;
		
		Submission(short depth, double s, unsigned int frq)
		{
			_depth = depth;
			_s = s;
			_freq = frq;
			_xn = 0; _yn = 0; _nn = 0;
		}
		Submission(vector<bool> * npath, short depth, double s, unsigned int frq)
		{
			_depth = depth;
			_s = s;
			_freq = frq;
			WedgeMatrix::FindIndices(npath, &_xn, &_yn, &_nn);
			_xn = depth - (_xn + 1);
		}
		Submission()
		{
			_depth = 0;
			_s = 0;
			_freq = 0;
			_xn = 0; _yn = 0; _nn = 0;
		}
	};
	
public:
	ThreadManager(int threads);
	virtual ~ThreadManager();
	
	void Spawn(short tid);
	void Kill(short tid);
	bool SafeToSpawn();
	void SetZWedge(WedgeMatrix * zw);
	short GetNextTID();
	
	void LongSubmit(short depth, vector<bool> * npath, float s, unsigned int freq, short tid);
	void ShortSubmit(short depth, float s, unsigned int freq, short tid);
	
protected:
	int _threads;
	int _activeThreads;
	WedgeMatrix * _zWedge;
	vector<Submission> ** _subs;
	bool * _tActive;
	
protected:
	void FeedToWedge(short tid);
	
public:
	static int SUB_THRESHOLD;
};

#endif // THREADMANAGER_HPP
