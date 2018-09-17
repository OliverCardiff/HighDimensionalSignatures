#include "Node.hpp"
#include <iostream>
#include <algorithm>
#include <math.h>
#include "Tree.hpp"
#include <thread>

using namespace std;

Node::Node(bool init)
{
	if(!init)
	{
		_scores.insert(pair<unsigned int, unsigned int>(0,0));
	}
}

Node::Node(unsigned int freq, bool init)
{
	_scores.insert(pair<unsigned int, unsigned int>(0,freq));
	
	/*if(init)
	{
		_children = new Node*[SETSIZE];
		for (unsigned int i = 0; i < SETSIZE; i++)
		{
			_children.at(i) = 0;
		}
	}*/
}

Node::~Node()
{
	_scores.clear();
	
	for (unsigned int i = 0; i < SETSIZE; i++)
	{
		if(_children.find(i) != _children.end())
		{
			delete _children.at(i);
		}
	}
	
	_children.clear();
}

void Node::KidsToTerminus(unsigned int ksID)
{
	unsigned int tks = ksID;
	for(short i = 0; i < (short)SETSIZE; i++)
	{
		if(i == 0)
		{
			tks = ksID * 2;
		}
		else
		{
			tks = ksID;
		}
		if(_children.find(i) != _children.end())
		{
			if(_children.at(i)->HasValue(tks))
			{
				unsigned int frq = _children.at(i)->KSFrq(tks);
				if(frq > 0)
				{
					TERMINUS[frq-1]++;
				}
			}
		}
	}
}

double Node::TerminalStruct()
{
	double agg = 0;
	
	for(double i = 0; i < 10000; i++)
	{
		if((i+1) > PLOIDY)
		{
			agg += (double)TERMINUS[i] * ceil(((i+1)/(double)PLOIDY)) * DEPTH;
		}
	}
	
	return agg;
}

double Node::PloidyRate()
{
	unsigned long int agg = 0;
	unsigned long int tracker = 0;
	
	for(int i = 0; i < 10000; i++)
	{
		agg += TERMINUS[i] * (i+1);
		tracker += TERMINUS[i];
	}
	
	unsigned long int pld = TERMINUS[PLOIDY - 1] * (PLOIDY - 1);
	unsigned long int paired = agg - tracker;
	
	if(paired > 0)
	{
		return ((double)pld / (double)paired) * 100;
	}
	else
	{
		return 0;
	}
}

double Node::PairedSeqsRate()
{
	unsigned long int agg = 0;
	unsigned long int tracker = 0;
	
	for(int i = 0; i < 10000; i++)
	{
		agg += TERMINUS[i] * (i+1);
		tracker += TERMINUS[i];
	}
	
	unsigned long int paired = agg - tracker;
	
	return (double)paired / (double)agg;
}

double Node::UnexplodedOrdinance()
{
	unsigned long int agg = 0;
	unsigned long int tracker = 0;
	
	for(int i = 0; i < 10000; i++)
	{
		agg += TERMINUS[i] * (i+1);
		tracker += TERMINUS[i];
	}
	
	unsigned long int start = agg;
	unsigned long int limit = agg * agg;
	unsigned long int actual = 0;
	
	for(int i = 0; i < 10000; i++)
	{
		actual += tracker * tracker;
		
		tracker -= TERMINUS[i];
	}
	
	return 1 - ((double)(actual - start) / (double)(limit - start));
	
}

bool Node::CanBeDeleted()
{
	if(_scores.size() == 1 && _scores.at(0) == 0)
	{
		return true;
	}
	return false;
}

void Node::InstantiateChildren()
{
	for(short i = 0; i < (short)SETSIZE; i++)
	{
		_children.insert(pair<short, Node*>(i,new Node(0, true)));
	}
}

void Node::RunZDFS(WedgeMatrix * zWedge)
{
	vector<bool> npath;
	
	zDFS(zWedge, 0, &npath, 0, 0);
	
	zWedge->ResolveParetos();
}

void Node::ClearAggMemory(unsigned int ksID)
{
	unsigned int tks = 0;
	for(short i = 0; i < (short)SETSIZE; i++)
	{
		if(i == 0)
		{
			tks = ksID * 2;
		}
		else
		{
			tks = ksID;
		}
		if(_children.find(i) != _children.end())
		{
			if(_children.at(i)->HasValue(tks))
			{
				_children.at(i)->ClearAggMemory(tks);
			}
			if(_children.at(i)->CanBeDeleted())
			{
				delete _children.at(i);
				_children.erase(i);
				//CHILDREN_KILLED++;
			}
		}
	}
	
	//int sz1 = _scores.size();
	_scores.erase(ksID);
	//int sz2 = _scores.size();
	//ID_DELETE += sz1 - sz2;
}

//retired
void Node::zDFS(WedgeMatrix * zWedge, int depth, vector<bool> * npath, double parentDistinction, unsigned int ksID)
{
	if(depth < (short)DEPTH)
	{
		if(Expanded(ksID) && depth > 0)
		{
			npath->push_back(true);
			_children.at(0)->zDFS(zWedge, depth+1, npath, parentDistinction, (ksID * 2) + 1);
		}
		
		double bd = FindBaseDistinction(ksID);
		
		double ad = parentDistinction * bd;
		
		parentDistinction = 1 - bd;
		
		if(ksID != 0)
		{
			zWedge->SendParetoZ(depth, npath, ad, KSFrq(ksID));
			zWedge->SendStructZ(depth, npath, ad, KSFrq(ksID));
		}
		else
		{
			zWedge->SendPurePareto(depth, ad, KSFrq(ksID));
			zWedge->SendPureStruct(depth, ad, KSFrq(ksID));
		}
		
		unsigned int tks = ksID;
		for(short i = 0; i < (short)SETSIZE; i++)
		{
			if(i == 0)
			{
				tks = ksID * 2;
			}
			else
			{
				tks = ksID;
			}
			if(_children.find(i) != _children.end())
			{
				if(_children.at(i)->HasValue(tks))
				{
					npath->push_back(false);
					_children.at(i)->zDFS(zWedge, depth+1, npath, parentDistinction, tks);
				}
			}
		}
	}
	npath->pop_back();
}

void Node::SendString(const char * kmer, short depth, int cov)
{
	_scores.at(0) += cov;
	char test = kmer[depth];
	
	bool terminalChildren = (depth != (short)DEPTH - 1);

	if (depth < (int)DEPTH)
	{
		for (unsigned int i = 0; i < SETSIZE; i++)
		{
			if (test == (*CHARSET)[i])
			{
				if (_children.find(i) != _children.end())
				{
					_children.at(i)->SendString(kmer, ++depth, cov);
				}
				else
				{
					_children.insert(pair<short, Node*>(i, new Node(0, terminalChildren)));
					_children.at(i)->SendString(kmer, ++depth, cov);
				}
				break;
			}
		}
	}
}

short Node::CountKids(unsigned int ksID)
{
	short kids = 0;
	
	if(_children.find(0) != _children.end())
	{
		if(_children.at(0)->HasValue(ksID * 2))
		{
			kids++;
		}
	}
	for(short i = 1; i < (short)SETSIZE; i++)
	{
		if(_children.find(i) != _children.end())
		{
			if(_children.at(i)->HasValue(ksID * 2))
			{
				kids++;
			}
		}
	}
	return kids;
}

bool Node::ExpandN(unsigned int ksID)
{
	short kids = 0;
	
	if(_children.find(0) != _children.end())
	{
		if(_children.at(0)->KSFrq(ksID * 2) > 0)
		{
			kids++;
		}
	}
	for (unsigned int i = 1; i < SETSIZE; i++)
	{
		if(_children.find(i) != _children.end())
		{
			if(_children.at(i)->KSFrq(ksID))
			{
				kids++;
			}
		}
	}
	if(kids > 1)
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool Node::Expanded(unsigned int ksID)
{
	if(_children.find(0) != _children.end())
	{
		unsigned int nks = (ksID * 2) + 1;
		if(_children.at(0)->HasValue(nks))
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	else
	{
		return false;
	}
}

bool Node::HasValue(unsigned int ksID)
{
	if(ksID == 0)
	{
		return true;
	}
	else if(_scores.find(ksID) != _scores.end())
	{
		return true;
	}
	return false;
}

void Node::MakeSureChildExists(short i)
{
	if(_children.find(i) == _children.end())
	{
		_children.insert(pair<short, Node*>(i, new Node(0,true)));
	}
}

//This has now been merged with zDFS
void Node::FindNextExpand(short depth, vector<bool> * npath, unsigned int ksID, short nCnt, double parentDistinction, WedgeMatrix * zWedge, short tid)
{
	bool doExpand = false;
	bool passDepth = false;
	bool expanded = false;
	
	if(depth < (int)DEPTH-1)
	{
		passDepth = true;
	}
	
	if(nCnt < NMAX && !Expanded(ksID))
	{
		if(ExpandN(ksID))
		{
			doExpand = true;
		}
	}
	
	unsigned int tks = ksID;
	
	if (passDepth)
	{
		if(nCnt < NMAX && doExpand && depth > 0)
		{
			if(_children.find(0) == _children.end())
			{
				_children.insert(pair<short, Node*>(0, new Node(0, true)));
			}
			unsigned int nks = (ksID * 2) + 1;
			_children.at(0)->CreateID(nks);
			_children.at(0)->PropKSFrq(KSFrq(ksID), nks);
			vector<KSNode> nodes;
			SubtreeAgg(ksID, depth, &nodes);
			expanded = true;
		}
	}
	
	if(expanded)
	{
		npath->push_back(true);
		zWedge->AddNFreq(_children.at(0)->KSFrq((ksID * 2) + 1));
		_children.at(0)->FindNextExpand(depth+1, npath, (ksID * 2) + 1, nCnt+1, parentDistinction, zWedge, tid);
	}
	
	//ks != 0 because backpropsubtract does the 0-ID leaves
	if(depth == (short)DEPTH - 1 && ksID != 0)
	{
		KidsToTerminus(ksID);
	}
	
	//The bit where the readings are actually taken
	if(depth < (short)DEPTH)
	{
		double bd = FindBaseDistinction(ksID);
			
		double ad = parentDistinction * bd;
		
		parentDistinction = 1 - bd;
		
		if(ad > 0)
		{
			if(ksID != 0)
			{
				Tree::TMAN->LongSubmit(depth, npath, ad, KSFrq(ksID), tid);
			}
			else
			{
				Tree::TMAN->ShortSubmit(depth, ad, KSFrq(ksID), tid);
			}
		}
	}
	
	vector<thread> thrs;
	vector<vector<bool> *> paths;
	vector<short> tids;
	short assigned = 0;
	short kids = CountKids(ksID);
	
	for (unsigned int i = 0; i < SETSIZE; i++)
	{
		if(i == 0)
		{
			tks = ksID * 2;
		}
		else
		{
			tks = ksID;
		}
		
		if(_children.find(i) != _children.end())
		{
			if (_children.at(i)->HasValue(tks))
			{
				bool madeThread = false;
				
				if(assigned < kids - 1 && Tree::TMAN->SafeToSpawn() && depth < 5)
				{
					short tid2 = Tree::TMAN->GetNextTID();
					if(tid2 != -1)
					{
						Tree::TMAN->Spawn(tid2);
						tids.push_back(tid2);
						vector<bool> * np2 = new vector<bool>();
						(*np2) = (*npath);
						np2->push_back(false);
						paths.push_back(np2);
						
						thrs.push_back(thread(&Node::FindNextExpand, _children.at(i), depth+1, np2, tks, nCnt, parentDistinction, zWedge, tid2));
						madeThread = true;
						//cout << "Making Thread: " << tid2 << "\n";
						
						assigned++;
					}
				}
				
				if(!madeThread)
				{
					npath->push_back(false);
					_children.at(i)->FindNextExpand(depth+1, npath, tks, nCnt, parentDistinction, zWedge, tid);
				}
			}
		}
	}
	npath->pop_back();
	
	for(short i = 0; i < assigned; i++)
	{
		thrs[i].join();
		delete paths[i];
		Tree::TMAN->Kill(tids[i]);
		//cout << "Killing Thread: " << tids[i] << "\n";
	}
	
	if(ksID != 0)
	{
		ClearAggMemory(ksID);
	}
}

void Node::SubtreeAgg(unsigned int ksID, short depth, vector<KSNode> * _nodes)
{
	if(depth < (short)DEPTH)
	{
		short nsz = _nodes->size();
		if(nsz == 0)
		{
			unsigned int tks = ksID;
			for (unsigned int i = 0; i < SETSIZE; i++)
			{
				if(i == 0)
				{
					tks = ksID * 2;
				}
				else
				{
					tks = ksID;
				}
				if(_children.find(i) != _children.end())
				{
					if(_children.at(i)->HasValue(tks))
					{
						(*_nodes).push_back(KSNode(_children.at(0), (ksID * 2) + 1));
						_children.at(i)->SubtreeAgg(tks, depth+1, _nodes);
					}
				}
			}
		}
		else
		{
			//The bit where we actually stitch the trees together
			
			vector<unsigned int> vals;
			
			unsigned int tks = ksID;
			for (short i = 0; i < (short)SETSIZE; i++)
			{
				if(i == 0)
				{
					tks = ksID * 2;
				}
				else
				{
					tks = ksID;
				}
				if(_children.find(i) != _children.end())
				{
					if(_children.at(i)->HasValue(tks))
					{
						vals.push_back(_children.at(i)->KSFrq(tks));
					}
					else
					{
						vals.push_back(0);
					}
				}
				else
				{
					vals.push_back(0);
				}
			}
			//propagating frequency values
			(*_nodes)[nsz-1]._ref->PropagateValues(&vals, (*_nodes)[nsz-1]._ksID);
			//then continue recursion
			
			unsigned int nks = (*_nodes)[nsz-1]._ksID;
			
			for (unsigned int i = 0; i < SETSIZE; i++)
			{
				if(i == 0)
				{
					tks = ksID * 2;
					nks = (*_nodes)[nsz-1]._ksID * 2;
				}
				else
				{
					tks = ksID;
					nks = (*_nodes)[nsz-1]._ksID;
				}
				if(_children.find(i) != _children.end())
				{
					if(_children.at(i)->HasValue(tks))
					{
						(*_nodes)[nsz-1]._ref->MakeSureChildExists((short)i);
						
						(*_nodes).push_back(KSNode((*_nodes)[nsz-1]._ref->GetChild(i), nks));
						_children.at(i)->SubtreeAgg(tks, depth+1, _nodes);
					}
				}
			}
		}
	}
	else
	{
		int nsz = (int)(*_nodes).size();
		unsigned int fq = (*_nodes)[nsz-1]._ref->KSFrq((*_nodes)[nsz-1]._ksID);
		if(fq > 1)
		{
			Node::TERMINUS[fq - 1]++;
			unsigned int penalty = fq - 1;
			
			for(int i = 0; i < nsz; i++)
			{
				(*_nodes)[i]._ref->SubKSFrq(penalty, (*_nodes)[i]._ksID);
			}
		}
	}
	_nodes->pop_back();
}

unsigned int Node::BackPropSubtract(short depth)
{
	unsigned int penalty = 0;
	
	if(depth == (short)DEPTH - 1)
	{
		KidsToTerminus(0);
	}
	
	if(depth < (short)DEPTH)
	{
		for (unsigned int i = 0; i < SETSIZE; i++)
		{
			if(_children.find(i) != _children.end())
			{
				penalty += _children.at(i)->BackPropSubtract(depth + 1);
			}
		}
	}
	else
	{
		unsigned int fq = KSFrq(0);
		if(fq > 1)
		{
			penalty += fq - 1;
		}
	}
	
	SetFreq(KSFrq(0) - penalty);
	
	return penalty;
}

//A Function to access the frequency of a given 
//aggregate/within node score, of a given annotation type
unsigned int Node::KSFrq(unsigned int ksID)
{	
	if(_scores.find(ksID) != _scores.end())
	{
		return _scores.at(ksID);
	}
	else
	{
		return 0;
	}
}

void Node::PropKSFrq(unsigned int freq, unsigned int ksID)
{
	if(_scores.find(ksID) == _scores.end())
	{
		CreateID(ksID);
		_scores.at(ksID) = freq;
	}
	else
	{
		_scores.at(ksID) += freq;
	}
}

void Node::SubKSFrq(unsigned int freq, unsigned int ksID)
{
	_scores.at(ksID) -= freq;
}

double Node::FindBaseDistinction(unsigned int ksID)
{
	double retVal = 0;
	int freq = KSFrq(ksID);
	
	if(freq > 1)
	{
		double fp = (double)freq;
		double fc = 0;
		double fm = ((freq/SETSIZE) * SETSIZE * SETSIZE) + ((freq % SETSIZE) * (freq % SETSIZE));
		double test = 0;
		
		//Finding Fc from here on out!
		vector<double> childF(SETSIZE, 0);
		
		if(_children.find(0) != _children.end())
		{	//Left turn carries all node doublings due to aggregation nodes doubling the K-space
			childF[0] = _children.at(0)->KSFrq(ksID*2);
			test += childF[0];
		}
		else
		{
			childF[0] = 0;
		}
		
		if(test == 0)
		{
			return 0;
		}
		
		for (unsigned short i = 1; i < SETSIZE; i++)
		{
			if(_children.find(i) != _children.end())
			{
				childF[i] = _children.at(i)->KSFrq(ksID);
				test += childF[i];
			}
			else
			{
				childF[i] = 0;
			}
		}
		
		sort(childF.begin(), childF.end());
		
		double accu = 0;
		
		for (unsigned short i = 1; i < SETSIZE; i++)
		{
			accu += (childF[i] - childF[i-1]) * (SETSIZE-i) * (SETSIZE-i);
		}
		
		fc = accu + (childF[0] * SETSIZE * SETSIZE);
		
		retVal = (fc - fp) / (fm - fp);
	}
	
	return retVal;
}

void Node::Increment(unsigned int ksID)
{	
	_scores.at(ksID)++;
}

void Node::SetFreq(unsigned int freq)
{
	_scores.at(0) = freq;
}

Node * Node::GetChild(short index)
{
	return _children.at(index);
}

void Node::CreateID(unsigned int ksID)
{
	//ID_CREATE++;
	_scores.insert(pair<unsigned int, unsigned int>(ksID, 0));
}

unsigned int Node::GetFreq()
{
	return _scores.at(0);
}

void Node::PropagateValues(vector<unsigned int> * vals, unsigned int ksID)
{
	unsigned int sz = vals->size();
	
	unsigned int tks = ksID;
	
	for(short i = 0; i < (short)sz; i++)
	{
		if(i == 0)
		{
			tks = ksID * 2;
		}
		else
		{
			tks = ksID;
		}
		if((*vals)[i])
		{
			if(_children.find(i) == _children.end())
			{
				_children.insert(pair<short, Node*>(i, new Node(0, true)));
			}
			_children.at(i)->PropKSFrq((*vals)[i], tks);
		}
	}
}

/* RETIRED
void Node::ExpandNMask(const char * kmer, short depth, unsigned int ksID, vector<unsigned int> * annoIDs, short currLeftTurn, short expands)
{
	bool doExpand = false;
	vector<unsigned int> expAnnots;
	bool terminal = false;
	bool increment = false;
	bool passDepth = false;
	
	if(depth < (int)DEPTH)
	{
		passDepth = true;
	}
	
	if(currLeftTurn == expands)
	{
		increment = true;
	}
	
	if(currLeftTurn >= NMAX)
	{
		terminal = true;
	}
	
	if(increment)
	{
		Increment(ksID, 0);
		
		for(unsigned int i = 0; i < annoIDs->size(); i++)
		{
			if(!HasValue(ksID, (*annoIDs)[i]))
			{
				_annoScores.at(ksID).insert(pair<unsigned short, unsigned int>((short)(*annoIDs)[i], 0));
			}
			
			Increment(ksID, (*annoIDs)[i]);
		}
	}
	else if(!terminal && passDepth)
	{
		if(!Expanded(ksID, 0))
		{
			if(ExpandN(ksID, 0))
			{
				doExpand = true;
				for(unsigned int i = 0; i < annoIDs->size(); i++)
				{
					if(!Expanded(ksID, (*annoIDs)[i]))
					{
						if(ExpandN(ksID, (*annoIDs)[i]))
						{
							expAnnots.push_back((*annoIDs)[i]);
						}
					}
				}
			}
		}
		else
		{
			for(unsigned int i = 0; i < annoIDs->size(); i++)
			{
				if(Expanded(ksID, (*annoIDs)[i]))
				{
					expAnnots.push_back((*annoIDs)[i]);
				}
			}
		}
	}
	else
	{
		return;
	}

	char test = kmer[depth];
	short tleft = currLeftTurn;
	short tks = ksID;
	
	if (passDepth)
	{
		if(doExpand)
		{
			if(_children.find(0) == _children.end())
			{
				_children.at(0) = new Node(0, true);
			}
			unsigned int nks = (ksID * 2) + 1;
			_children.at(0)->CreateID(nks);
			_children.at(0)->ExpandNMask(kmer, depth+1, nks, &expAnnots, tleft+1, expands);
		}
		else if(Expanded(ksID,0))
		{
			unsigned int nks = (ksID * 2) + 1;
			_children.at(0)->ExpandNMask(kmer, depth+1, nks, &expAnnots, tleft+1, expands);
		}
		for (unsigned int i = 0; i < SETSIZE; i++)
		{
			if (test == (*CHARSET)[i])
			{
				if(i == 0)
				{
					tleft = currLeftTurn;
					if(!terminal)
					{
						tks = ksID * 2;
					}
				}
				else
				{
					tleft = currLeftTurn;
					tks = ksID;
				}
				if(increment)
				{
					if(_children.find(i) == _children.end())
					{
						_children.at(i) = new Node(0, true);
					}
					if (!_children.at(i)->HasValue(tks, 0))
					{
						_children.at(i)->CreateID(tks);
					}
					_children.at(i)->ExpandNMask(kmer, depth+1, tks, annoIDs, tleft, expands);
				}
				else
				{
					_children.at(i)->ExpandNMask(kmer, depth+1, tks, annoIDs, tleft, expands);
				}
			}
		}
	}
}
*/