#include "WedgeMatrix.hpp"
#include "Node.hpp"
#include <math.h>

WedgeMatrix::WedgeMatrix(int k, bool full)
{
	_k = k;
	_full = full;
	sumF = 0; sumN = 0; sumP = 0; sumS = 0;
	_NFReqs = 0;
	_hasNull = false;
	_isNull = false;
	_headFreq = 0;
	_isLocked = false;
	
	if(full)
	{
		_sMemory = new vector<float>***[_k];
		_fMemory = new vector<int>***[_k];
		_zFMemory = new vector<float>*[_k];
		_zNMemory = new vector<int>*[_k];
	}

	_fTracker = new unsigned int**[_k];
	_dfTracker = new unsigned int**[_k];
	_nTracker = new unsigned int**[_k];
	
	_nMemoryF = new vector<float>*[_k];
	_nMemoryS = new vector<float>*[_k];
	_nMemoryShape = new vector<float>*[_k];
	_nMemoryPN = new vector<float>*[_k];
	
	_sMeans = new float**[_k];
	_sDeviations = new float**[_k];
	
	for(int i = 0; i < _k; i++)
	{
		_zStruct.push_back(0);
		_zDev.push_back(0);
		
		if(full)
		{
			_sMemory[i] = new vector<float>**[i+1];
			_fMemory[i] = new vector<int>**[i+1];
			_zFMemory[i] = new vector<float>();
			_zNMemory[i] = new vector<int>();
		}

		_zFTracker.push_back(0);
		_zNTracker.push_back(0);
		
		_nMemoryF[i] = new vector<float>(Node::NMAX+1, 0);
		_nMemoryS[i] = new vector<float>(Node::NMAX+1, 0);
		_nMemoryShape[i] = new vector<float>(Node::NMAX+1, 0);
		_nMemoryPN[i] = new vector<float>(Node::NMAX+1, 0);
		
		_fTracker[i] = new unsigned int*[i+1];
		_dfTracker[i] = new unsigned int*[i+1];
		_nTracker[i] = new unsigned int*[i+1];
		
		_sMeans[i] = new float*[i+1];
		_sDeviations[i] = new float*[i+1];
		
		for(int j = 0; j <= i; j++)
		{
			if(full)
			{
				_sMemory[i][j] = new vector<float>*[j+1];
				_fMemory[i][j] = new vector<int>*[j+1];
			}
			_fTracker[i][j] = new unsigned int[j+1];
			_dfTracker[i][j] = new unsigned int[j+1];
			_nTracker[i][j] = new unsigned int[j+1];
			
			_sMeans[i][j] = new float[j+1];
			_sDeviations[i][j] = new float[j+1];
			
			for(int l = 0; l <= j; l++)
			{
				if(full)
				{
					_sMemory[i][j][l] = new vector<float>();
					_fMemory[i][j][l] = new vector<int>();
				}

				_fTracker[i][j][l] = 0;
				_dfTracker[i][j][l] = 0;
				_nTracker[i][j][l] = 0;

				_sMeans[i][j][l] = 0;
				_sDeviations[i][j][l] = 0;
			}
		}
	}
}

WedgeMatrix::~WedgeMatrix()
{
	for(int i = 0; i < _k; i++)
	{		
		for(int j = 0; j <= i; j++)
		{			
			delete[] _sMeans[i][j];
			delete[] _sDeviations[i][j];

			delete[] _fTracker[i][j];
			delete[] _nTracker[i][j];

		}

		delete[] _fTracker[i];
		delete[] _nTracker[i];

		delete[] _sMeans[i];
		delete[] _sDeviations[i];
	}

	delete[] _fTracker;
	delete[] _nTracker;

	delete[] _sMeans;
	delete[] _sDeviations;
}

void WedgeMatrix::Lock()
{
	_isLocked = true;
}

void WedgeMatrix::Unlock()
{
	_isLocked = false;
}

bool WedgeMatrix::IsLocked()
{
	return _isLocked;
}

void WedgeMatrix::MergeIn()
{
	for(int i = 0; i < _k; i++)
	{
		_sMeans[i][0][0] = _zStruct[i];
		_sDeviations[i][0][0] = _zDev[i];
		_fTracker[i][0][0] = _zFTracker[i];
		_nTracker[i][0][0] = _zNTracker[i];
	}
}

unsigned int WedgeMatrix::GetFreqM(short i, short j, short k)
{
	return _fTracker[i][j][k];
}

unsigned int WedgeMatrix::GetShapeN(short i, short j, short k)
{
	return _nTracker[i][j][k];
}

unsigned int WedgeMatrix::GetDepthFreqM(short i, short j, short k)
{
	return _dfTracker[i][j][k];
}

void WedgeMatrix::AddNFreq(unsigned int f)
{
	_NFReqs += f;
}

void WedgeMatrix::AddHeadFreq(unsigned int f)
{
	_headFreq = f;
}

void WedgeMatrix::SetAsNull()
{
	_isNull = true;
}

bool WedgeMatrix::IsNull()
{
	return _isNull;
}

double WedgeMatrix::GetNormalisedS()
{
	return sumS / (((double)_headFreq * (double)_k) + (double)_NFReqs * (double)_k);
}

double WedgeMatrix::GetSumP()
{
	return sumN / sumP;
}

double WedgeMatrix::GetSumStruct()
{
	return sumS;
}

double WedgeMatrix::GetSumFs()
{
	return (((double)_headFreq * (double)_k) + (double)_NFReqs * (double)_k);
}

double WedgeMatrix::GetParetoN()
{
	return sumN;
}

float WedgeMatrix::GetZScoreD(short i)
{
	return _zDev[i];
}

float WedgeMatrix::GetZScoreM(short i)
{
	return _zStruct[i];
}

double WedgeMatrix::CorrectFreq(unsigned int frq, short depth)
{
	return (double)frq / max(1.0, (double)_headFreq/(double)pow(SETSIZE, depth));
}

void WedgeMatrix::NormaliseStructs()
{
	for(int i = 0; i < _k; i++)
	{
		if(_zFTracker[i] > 0)
		{
			_zStruct[i] = _zStruct[i] / _zFTracker[i];
		}
		
		for(int j = 0; j <= i; j++)
		{
			for(int l = 0; l <= j; l++)
			{
				if(_fTracker[i][j][l] > 0)
				{
					_sMeans[i][j][l] = _sMeans[i][j][l] / _fTracker[i][j][l];
				}
			}
		}
	}
}

void WedgeMatrix::ReWeightStructs()
{
	for(int i = 0; i < _k; i++)
	{
		if(_zFTracker[i] > 0)
		{
			_zStruct[i] = _zStruct[i] * _zFTracker[i];
		}
		
		for(int j = 0; j <= i; j++)
		{
			for(int l = 0; l <= j; l++)
			{
				if(_fTracker[i][j][l] > 0)
				{
					_sMeans[i][j][l] = _sMeans[i][j][l] * _fTracker[i][j][l];
				}
			}
		}
	}
}

float WedgeMatrix::GetScoreD(short i, short j, short k)
{
	return _sDeviations[i][j][k];
}

float WedgeMatrix::GetScoreM(short i, short j, short k)
{
	return _sMeans[i][j][k];
}

float WedgeMatrix::GetNScoreS(short d, short n)
{
	return (*_nMemoryS[d])[n];
}

float WedgeMatrix::GetNScoreF(short d, short n)
{
	return (*_nMemoryF[d])[n];
}

float WedgeMatrix::GetNScoreShape(short d, short n)
{
	return (*_nMemoryShape[d])[n];
}

float WedgeMatrix::GetNScorePN(short d, short n)
{
	return (*_nMemoryPN[d])[n];
}

void WedgeMatrix::SetNull(WedgeMatrix * nullMat)
{
	_referenceNull = nullMat;
	if(nullMat)
	{
		_hasNull = true;
	}
	else
	{
		_hasNull = false;
	}
}

void WedgeMatrix::MinZero(float * val)
{
	if((*val) < 0)
	{
		(*val) = 0;
	}
}

void WedgeMatrix::SendPureStruct(short depth, float s, unsigned int freq)
{
	double f2 = CorrectFreq(freq, depth);
	float p = s * f2 * depth;
	
	//if(p > PLOIDY)
	//{
		_zStruct[depth] += p;
		_zFTracker[depth] += f2;
		_dfTracker[depth][0][0] += f2 * depth;
		
		sumS += p;
		sumF += f2;
	//}
}

void WedgeMatrix::SendPurePareto(short depth, float s, unsigned int freq)
{
	double f2 = CorrectFreq(freq, depth);
	float p = s * f2;// * depth;
	
	if(s * f2 > 1) //PLOIDY)
	{
		_zNTracker[depth]++;
		_zDev[depth] += log(p);// - log(PLOIDY);
		
		sumN += 1;
		sumP += log(p);// - log(PLOIDY);
	}
}

void WedgeMatrix::ResolveParetos()
{
	if(!_hasNull)
	{
		sumP = 0;
		sumN = 0;
		for(int i = 0; i < _k; i++)
		{
			for(int j = 0; j <= i; j++)
			{
				for(int l = 0; l <= j; l++)
				{
					if(_sDeviations[i][j][l] > 0)
					{
						sumP += _sDeviations[i][j][l];
						sumN += _nTracker[i][j][l];
						//_sDeviations[i][j][l] = _nTracker[i][j][l] / _sDeviations[i][j][l];
					}
				}
			}
		}
	}
	else
	{
		sumP = 0;
		sumN = 0;
		for(int i = 0; i < _k; i++)
		{
			for(int j2 = 0; j2 <= Node::NMAX; j2++)
			{
				float zref = _referenceNull->GetNScoreS(i, j2);
				
				if((*_nMemoryS[i])[j2] > 0)
				{
					(*_nMemoryShape[i])[j2] += (*_nMemoryPN[i])[j2] * log((*_nMemoryS[i])[j2]/((*_nMemoryS[i])[j2] + zref));
					if((*_nMemoryShape[i])[j2] <= 0)
					{
						(*_nMemoryPN[i])[j2] = max(0.0, (*_nMemoryPN[i])[j2] - ceil(((*_nMemoryShape[i])[j2] * -1) 
							/ ((log((*_nMemoryS[i])[j2] / ((*_nMemoryS[i])[j2] + zref)) * -1))));
						(*_nMemoryShape[i])[j2] = 0;
					}
				}
				else
				{
					(*_nMemoryShape[i])[j2] = 0;
					(*_nMemoryPN[i])[j2] = 0;
				}
			}
			for(int j = 0; j <= i; j++)
			{
				for(int l = 0; l <= j; l++)
				{
					if(_sDeviations[i][j][l] > 0)
					{
						float zref = _referenceNull->GetScoreM(i, j, l);
						if(_sMeans[i][j][l] > 0)
						{
							_sDeviations[i][j][l] += _nTracker[i][j][l] * log(_sMeans[i][j][l]/(_sMeans[i][j][l] + zref));
							if(_sDeviations[i][j][l] > 0)
							{
								sumP += _sDeviations[i][j][l];
								sumN += _nTracker[i][j][l];
								//_sDeviations[i][j][l] = _nTracker[i][j][l] / _sDeviations[i][j][l];
							}
							else
							{
								//CALCULATE MINIMUM N REDUCTION TO GET BACK TO ZERO
								_nTracker[i][j][l] = max(0.0, _nTracker[i][j][l] - ceil((_sDeviations[i][j][l] * -1) 
								/ ((log(_sMeans[i][j][l] / (_sMeans[i][j][l] + zref)) * -1))));
								sumN += _nTracker[i][j][l];
								_sDeviations[i][j][l] = 0;
							}
						}
						else
						{
							_sDeviations[i][j][l] = 0;
							_nTracker[i][j][l] = 0;
						}
					}
				}
			}
		}
	}
}

void WedgeMatrix::FindIndices(vector<bool> * npath, short * xn, short * yn, short * nn)
{
	(*xn) = -1;
	(*yn) = -1;
	(*nn) = 0;
	
	for(unsigned short i = 1; i < npath->size(); i++)
	{
		if((*npath)[i])
		{
			(*nn)++;
			if((*xn) == -1)
			{
				(*xn) = i-1;
			}
			(*yn) = i-1;
		}
	}
	
	(*yn) -= (*xn);
}

void WedgeMatrix::SendStructZ(short depth, vector<bool> * npath, float s, unsigned int freq)
{
	short xn; short yn; short nn;
	FindIndices(npath, &xn, &yn, &nn);
	xn = depth - (xn + 1);
	double f2 = CorrectFreq(freq, (depth - nn));
	float p = s * f2 * (depth - nn);
	
	_sMeans[depth][xn][yn] += p;
	_fTracker[depth][xn][yn] += f2;
	_dfTracker[depth][xn][yn] += f2 * (depth - nn);
	
	sumS += p;
	sumF += f2;
}

void WedgeMatrix::SendParetoZ(short depth, vector<bool> * npath, float s, unsigned int freq)
{
	short xn; short yn; ; short nn;
	FindIndices(npath, &xn, &yn, &nn);
	xn = depth - (xn + 1);
	float p = s;// * (depth - nn);
	double f2 = CorrectFreq(freq, (depth - nn));
	
	if(s * f2 > 1) //PLOIDY)
	{
		_sDeviations[depth][xn][yn] += log(p * f2);// - log(PLOIDY);
		_nTracker[depth][xn][yn]++;
		
		sumN += 1;
		sumP += log(p * f2);
	}
}

void WedgeMatrix::TSendStructZ(short depth, short xn, short yn, short nn, float s, unsigned int freq)
{
	double f2 = CorrectFreq(freq, (depth - nn));
	float p = s * f2 * (depth - nn);
	
	_sMeans[depth][xn][yn] += p;
	_fTracker[depth][xn][yn] += f2;
	_dfTracker[depth][xn][yn] += f2 * (depth - nn);
	
	sumS += p;
	sumF += f2;
}

void WedgeMatrix::TSendParetoZ(short depth, short xn, short yn, short nn, float s, unsigned int freq)
{
	float p = s;// * (depth - nn);
	double f2 = CorrectFreq(freq, (depth - nn));
	
	if(s * f2 > 1) //PLOIDY)
	{
		_sDeviations[depth][xn][yn] += log(p * f2);// - log(PLOIDY);
		_nTracker[depth][xn][yn]++;
		
		sumN += 1;
		sumP += log(p * f2);
	}
}

void WedgeMatrix::TSendBoth(short depth, short xn, short yn, short nn, float s, unsigned int freq)
{
	float p = s;// * (depth - nn);
	double f2 = CorrectFreq(freq, (depth - nn));
	
	if(s * f2 > 1) //PLOIDY)
	{
		_sDeviations[depth][xn][yn] += log(p * f2);// - log(PLOIDY);
		_nTracker[depth][xn][yn]++;
		
		(*_nMemoryShape[depth])[nn] += log(p * f2);
		(*_nMemoryPN[depth])[nn]++;
		
		sumN += 1;
		sumP += log(p * f2);
	}
	
	p = s * f2 * (depth - nn);
	
	_sMeans[depth][xn][yn] += p;
	_fTracker[depth][xn][yn] += f2;
	_dfTracker[depth][xn][yn] += f2 * (depth - nn);
	
	(*_nMemoryS[depth])[nn] += s * f2;
	(*_nMemoryF[depth])[nn] += f2;
	
	sumS += p;
	sumF += f2;
	
	if(_full)
	{
		if(f2 >= 1)
		{
			(*_sMemory[depth][xn][yn]).push_back(s);
			(*_fMemory[depth][xn][yn]).push_back((int)f2);
		}
	}
	
}

void WedgeMatrix::DeleteMemory()
{
	for(int i = 0; i < _k; i++)
	{		
		for(int j = 0; j <= i; j++)
		{			
			for(int l = 0; l <= j; l++)
			{
				delete _sMemory[i][j][l];
				delete _fMemory[i][j][l];
				
			}
			delete[] _sMemory[i][j];
			delete[] _fMemory[i][j];
		}
		delete[] _sMemory[i];
		delete[] _fMemory[i];
		delete _zFMemory[i];
		delete _zNMemory[i];
	}
	delete[] _sMemory;
	delete[] _fMemory;
	delete[] _zFMemory;
	delete[] _zNMemory;
}

void WedgeMatrix::SubtractNull(bool doSD)
{
	if(_hasNull)
	{
		sumS -= _referenceNull->GetSumStruct();
		if(sumS < 0) sumS = 0;
		for(int i = 0; i < _k; i++)
		{
			_zStruct[i] -= _referenceNull->GetZScoreM(i);
			MinZero(&_zStruct[i]);
			
			for(int j2 = 0; j2 <= Node::NMAX; j2++)
			{
				(*_nMemoryS[i])[j2] -= _referenceNull->GetNScoreS(i, j2);
				
				float refFreq = _referenceNull->GetNScoreF(i, j2);
				if(refFreq < (*_nMemoryF[i])[j2])
				{
					(*_nMemoryF[i])[j2] -= refFreq;
				}
				else
				{
					(*_nMemoryF[i])[j2] = 0;
				}
				
				MinZero(&(*_nMemoryS[i])[j2]);
				MinZero(&(*_nMemoryF[i])[j2]);
			}
			
			if(doSD)
			{
				_zDev[i] -= _referenceNull->GetZScoreD(i);
				MinZero(&_zDev[i]);
			}
			for(int j = 0; j <= i; j++)
			{			
				for(int l = 0; l <= j; l++)
				{
					_sMeans[i][j][l] -= _referenceNull->GetScoreM(i,j,l);
					
					unsigned int refDF = _referenceNull->GetDepthFreqM(i,j,l);
					unsigned int refF = _referenceNull->GetFreqM(i,j,l);
					
					if(refDF < _dfTracker[i][j][l])
					{
						_dfTracker[i][j][l] -= refDF;
					}
					else
					{
						_dfTracker[i][j][l] = 0;
					}
					if(refF < _fTracker[i][j][l])
					{
						_fTracker[i][j][l] -= refF;
					}
					else
					{
						_fTracker[i][j][l] = 0;
					}
					
					MinZero(&_sMeans[i][j][l]);
					
					if(doSD)
					{
						_sDeviations[i][j][l] -= _referenceNull->GetScoreD(i,j,l);
						MinZero(&_sDeviations[i][j][l]);
					}
				}
			}
		}
	}
}

void WedgeMatrix::ResolveFull()
{
	//TODO: Add single vector zFMemory and zNMemory resolve too (same function different vectors!!)
	for(int i = 0; i < _k; i++)
	{
		FindMSD(_zFMemory[i], _zNMemory[i], &(_zStruct[i]), &(_zDev[i]));
		
		for(int j = 0; j <= i; j++)
		{			
			for(int l = 0; l <= j; l++)
			{
				FindMSD(_sMemory[i][j][l], _fMemory[i][j][l], &_sMeans[i][j][l], &_sDeviations[i][j][l]);
			}
		}
	}
	DeleteMemory();
}

void WedgeMatrix::FindMSD(vector<float> * smem, vector<int> * fmem, float * meanVal, float * sdVal)
{
	float saccu = 0;
	float faccu = 0;
	float sScore = 0;	
	float wSD = 0;

	int sz = smem->size();
	
	for(int m = 0; m < sz; m++)
	{
		saccu += (*smem)[m] * (float)(*fmem)[m];
		faccu += (float)(*fmem)[m];
	}
	
	if(sz != 0)
	{
		sScore = saccu/faccu;
		if(sz > 1)
		{
			saccu = 0;
			float diff = 0;
			for(int m = 0; m < sz; m++)
			{
				diff = (*smem)[m] - sScore;
				saccu += diff * diff * (float)(*fmem)[m];
			}
			
			wSD = saccu / (((sz - 1) * faccu) / sz);
		}
	}
	(*meanVal) = sScore;
	(*sdVal) = wSD;
}