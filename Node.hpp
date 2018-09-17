#ifndef NODE_H
#define NODE_H

#include <string>
#include <vector>
#include <map>
#include "WedgeMatrix.hpp"

using namespace std;

class Node
{
	struct KSNode
	{
		unsigned int _ksID;
		Node * _ref;
		KSNode()
		{
			_ksID = 0;
			_ref = 0;
		}
		KSNode(Node * ref, unsigned int Id)
		{
			_ksID = Id;
			_ref = ref;
		}
	};
public:
	Node(bool init);
	Node(unsigned int freq, bool init);
	~Node();

public:
	void RunZDFS(WedgeMatrix * zWedge);
	void zDFS(WedgeMatrix * zWedge, int depth, vector<bool> * npath, double parentDistinction, unsigned int ksID);
	double FindBaseDistinction(unsigned int ksID);
	unsigned int KSFrq(unsigned int ksID);
	void PropKSFrq(unsigned int freq, unsigned int ksID);
	void Increment(unsigned int ksID);
	
	void SendString(const char * kmer, short depth, int cov);
	bool ExpandN(unsigned int ksID);
	bool Expanded(unsigned int ksID);
	unsigned int BackPropSubtract(short depth);
	void FindNextExpand(short depth, vector<bool> * npath, unsigned int ksID, short nCnt, double parentDistinction, WedgeMatrix * zWedge, short tid);
	void SubtreeAgg(unsigned int ksID, short depth, vector<KSNode> * _nodes);
	bool HasValue(unsigned int ksID);
	
	short CountKids(unsigned int ksID);
	
public:
	void InstantiateChildren();
	Node * GetChild(short index);
	void SetFreq(unsigned int freq);
	void SubKSFrq(unsigned int freq, unsigned int ksID);
	unsigned int GetFreq();
	void CreateID(unsigned int ksID);
	void PropagateValues(vector<unsigned int> * vals, unsigned int ksID);
	void MakeSureChildExists(short i);
	void ClearAggMemory(unsigned int ksID);
	bool CanBeDeleted();
	void KidsToTerminus(unsigned int ksID);
	
public:
	static unsigned int PLOIDY;
	static unsigned int SETSIZE;
	static std::vector<char> * CHARSET;
	static unsigned int DEPTH;
	static unsigned int AGGS;
	static char * PREFIX;
	static unsigned short NMAX;
	
	static vector<unsigned int> TERMINUS;
	
	static double UnexplodedOrdinance();
	static double PairedSeqsRate();
	static double PloidyRate();
	static double TerminalStruct();
	
protected:
	map<short, Node*> _children;
	map<unsigned int, unsigned int> _scores;
};

#endif // TREE_H