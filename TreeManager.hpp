#ifndef TREEMANAGER_HPP
#define TREEMANAGER_HPP

#include "Tree.hpp"
#include <iostream>
#include <fstream>
#include "WedgeMatrix.hpp"

class TreeManager
{
public:
	TreeManager(short depth, bool subtract, bool dnull);
	virtual ~TreeManager();
	
public:
	void RunZDSCS(char * genome, bool hasDepths, char * depthFile);
	void SubtractAndResolve();
	void ZWedgeSummary();
	void CollapseSDs();
	void WriteSDMatrices(char * prefix);
	void WriteWedgeFiles(char * prefix, bool subtracted);
	
protected:
	short _depth;
	bool _doSubtract;
	bool _doNull;
	Tree * _mainTree;
	Tree * _nullTree;
	WedgeMatrix * _nullWedge;
	WedgeMatrix * _mainWedge;
};

#endif // TREEMANAGER_HPP
