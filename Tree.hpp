#ifndef TREE_H
#define TREE_H
#include "Node.hpp"
#include "WedgeMatrix.hpp"
#include <string>
#include <vector>
#include <algorithm>
#include "ThreadManager.hpp"

class Tree
{
public:
	Tree(int depth, bool isNull, bool dosubs);
	~Tree();
	void ReadFasta(char * genome);
	void ReadDepths(char * depthFile);
	void ExpandTree();
	void SetSkip(double skip);
	void FirstTreeSummary();
	
	double GetHeadFreq();
	double GetSEscape();
	void GetZeroS();
	void StateZSummary(bool escapes);
	void ResolveParetos();
	WedgeMatrix * GetZWedge();
	void SetNullReference(WedgeMatrix * nullMat);
	
public:
	//Fingerprint class
	static ThreadManager * TMAN;
	
protected:
	Node * _root;
	bool _subtract;
	bool _hasDepths;
	int _depth;
	void ReadInString(std::string * strs, int cov);
	void ReadNullString(string * strs, int cov);
	void RevComp(std::string * strs);
	bool _null;
	WedgeMatrix * _zWedge;
	vector<int> _parallelDepths;
	double _ordinance;
	double _psRate;
	double _plRate;
	double _escapeS;
	double _skipRate;
	unsigned int _maxPenalty;
};

#endif // TREE_H
