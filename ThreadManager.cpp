#include "ThreadManager.hpp"

ThreadManager::ThreadManager(int threads)
{
	_threads = threads;
	_activeThreads = 1;
	_subs = new vector<Submission>*[_threads];
	_tActive = new bool[_threads];
	_tActive[0] = true;
	
	for(int i = 0; i < threads; i++)
	{
		if(i > 0)
		{
			_tActive[i] = false;
		}
		_subs[i] = new vector<Submission>(SUB_THRESHOLD);
	}
}

ThreadManager::~ThreadManager()
{
	for(int i = 0; i < _threads; i++)
	{
		delete _subs[i];
	}
	delete [] _subs;
}

void ThreadManager::Spawn(short tid)
{
	_tActive[tid] = true;
	_activeThreads++;
}

void ThreadManager::Kill(short tid)
{
	_tActive[tid] = false;
	_activeThreads--;
}

short ThreadManager::GetNextTID()
{
	for(short i = 0; i < _threads; i++)
	{
		if(!_tActive[i])
		{
			return i;
		}
	}
	return -1;
}

void ThreadManager::SetZWedge(WedgeMatrix * zw)
{
	_zWedge = zw;
}

bool ThreadManager::SafeToSpawn()
{
	if(_activeThreads < _threads)
	{
		return true;
	}
	return false;
}

void ThreadManager::ShortSubmit(short depth, float s, unsigned int freq, short tid)
{
	_subs[tid]->push_back(Submission(depth, s, freq));
	
	FeedToWedge(tid);
}

void ThreadManager::LongSubmit(short depth, vector<bool> * npath, float s, unsigned int freq, short tid)
{
	_subs[tid]->push_back(Submission(npath, depth, s, freq));
	
	FeedToWedge(tid);
}

void ThreadManager::FeedToWedge(short tid)
{
	int slen = _subs[tid]->size();
	
	if(slen > SUB_THRESHOLD && !_zWedge->IsLocked())
	{
		_zWedge->Lock();
		for(int i = 0; i < slen; i++)
		{
			Submission * sb = &(*_subs[tid])[i];
			
			//_zWedge->TSendParetoZ(sb->_depth, sb->_xn, sb->_yn, sb->_nn, sb->_s, sb->_freq);
			//_zWedge->TSendStructZ(sb->_depth, sb->_xn, sb->_yn, sb->_nn, sb->_s, sb->_freq);
			_zWedge->TSendBoth(sb->_depth, sb->_xn, sb->_yn, sb->_nn, sb->_s, sb->_freq);
		}
		_subs[tid]->clear();
		_zWedge->Unlock();
	}
}