#include "TreeManager.hpp"
#include "Node.hpp"

using namespace std;

TreeManager::TreeManager(short depth, bool subtract, bool dnull)
{
	_doNull = dnull;
	_depth = depth;
	_doSubtract = subtract;
}

TreeManager::~TreeManager()
{
}

void TreeManager::ZWedgeSummary()
{
	cout <<  "\n~ ~ ~ Z-DSCS Wedge Summary ~ ~ ~\n";
	
	for(int i = 0; i < _depth; i++)
	{
		cout << "i: " << i;
		for(int j = 0; j <= i; j++)
		{			
			cout << " j: " << j;
			for(int l = 0; l <= j; l++)
			{
				cout << "\t" << _mainWedge->GetScoreM(i,j,l);
			}
		}
		cout << "\n\n";
	}
}

void TreeManager::CollapseSDs()
{
	_mainWedge->ResolveFull();
}

void TreeManager::WriteSDMatrices(char * prefix)
{
	string pref(prefix);
	string mainOut = pref + "_main_SD.txt";
	
	ofstream ofs_main(mainOut.c_str(), ofstream::out);
	
	ofs_main << "i\tj\tk\tStruct\tFreq\tMSD\n";
	
	for(int i = 0; i < _depth; i++)
	{
		if(i == 0)
		{
			ofs_main << i << "\t" << 0 << "\t" << 0 << "\t" 
			<< _mainWedge->GetScoreM(i, 0, 0) << "\t"
			<< _mainWedge->GetFreqM(i, 0, 0) << "\t"
			<< _mainWedge->GetScoreD(i, 0, 0) << "\n";
		}
		
		for(int j = 0; j < i; j++)
		{
			if(j == 0)
			{
				ofs_main << i << "\t" << j << "\t" << 0 << "\t" 
				<< _mainWedge->GetScoreM(i, j, 0) << "\t"
				<< _mainWedge->GetFreqM(i, j, 0) << "\t"
				<< _mainWedge->GetScoreD(i, j, 0) << "\n";
			}
			for(int l = 0; l < j; l++)
			{
				ofs_main << i << "\t" << j << "\t" << l << "\t" 
				<< _mainWedge->GetScoreM(i, j, l) << "\t"
				<< _mainWedge->GetFreqM(i, j, l) << "\t"
				<< _mainWedge->GetScoreD(i, j, l) << "\n";
			}
		}
	}
	ofs_main.close();
}

void TreeManager::WriteWedgeFiles(char * prefix, bool subtracted)
{
	string pref(prefix);
	string mainOut = pref + "_main";
	string nullOut = pref + "_null.txt";
	string mainCompact = pref + "_short_main";
	string nullCompact = pref + "_short_null.txt";
	
	if(subtracted)
	{
		mainOut = mainOut + "_sub.txt";
		mainCompact = mainCompact + "_sub.txt";
	}
	else
	{
		mainOut = mainOut + "_prior.txt";
		mainCompact = mainCompact + "_prior.txt";
	}
	
	ofstream ofs_mshrt(mainCompact.c_str(), ofstream::out);
	
	ofs_mshrt << "Depth\tNs\tStruct\tFreq\tShape\tShapeN\n";
	
	for(int i = 0; i < _depth; i++)
	{
		for(short j2 = 0; j2 <= (short)Node::NMAX; j2++)
		{
			ofs_mshrt << i << "\t" << j2 << "\t" 
			<< _mainWedge->GetNScoreS(i, j2) << "\t"
			<< _mainWedge->GetNScoreF(i, j2) << "\t"
			<< _mainWedge->GetNScoreShape(i, j2) << "\t"
			<< _mainWedge->GetNScorePN(i, j2) << "\n";
		}
	}
	
	ofs_mshrt.close();
	
	ofstream ofs_main(mainOut.c_str(), ofstream::out);
	
	ofs_main << "i\tj\tk\tStruct\tFreq\tDepthF\tShape\tN\n";
	
	for(int i = 0; i < _depth; i++)
	{
		if(i == 0)
		{
			ofs_main << i << "\t" << 0 << "\t" << 0 << "\t" 
			<< _mainWedge->GetScoreM(i, 0, 0) << "\t"
			<< _mainWedge->GetFreqM(i, 0, 0) << "\t"
			<< _mainWedge->GetDepthFreqM(i, 0, 0) << "\t"
			<< _mainWedge->GetScoreD(i, 0, 0) << "\t" 
			<< _mainWedge->GetShapeN(i, 0, 0) << "\n";
		}
		
		for(int j = 0; j < i; j++)
		{
			if(j == 0)
			{
				ofs_main << i << "\t" << j << "\t" << 0 << "\t" 
				<< _mainWedge->GetScoreM(i, j, 0) << "\t"
				<< _mainWedge->GetFreqM(i, j, 0) << "\t"
				<< _mainWedge->GetDepthFreqM(i, j, 0) << "\t"
				<< _mainWedge->GetScoreD(i, j, 0) << "\t" 
				<< _mainWedge->GetShapeN(i, j, 0) << "\n";
			}
			for(int l = 0; l < j; l++)
			{
				ofs_main << i << "\t" << j << "\t" << l << "\t" 
				<< _mainWedge->GetScoreM(i, j, l) << "\t"
				<< _mainWedge->GetFreqM(i, j, l) << "\t"
				<< _mainWedge->GetDepthFreqM(i, j, l) << "\t"
				<< _mainWedge->GetScoreD(i, j, l) << "\t" 
				<< _mainWedge->GetShapeN(i, j, l) << "\n";
			}
		}
	}
	ofs_main.close();
	
	if(subtracted)
	{
		ofstream ofs_nshrt(nullCompact.c_str(), ofstream::out);
	
		ofs_nshrt << "Depth\tNs\tStruct\tFreq\tShape\tShapeN\n";
		
		for(int i = 0; i < _depth; i++)
		{
			for(short j2 = 0; j2 <= (short)Node::NMAX; j2++)
			{
				ofs_nshrt << i << "\t" << j2 << "\t" 
				<< _nullWedge->GetNScoreS(i, j2) << "\t"
				<< _nullWedge->GetNScoreF(i, j2) << "\t"
				<< _nullWedge->GetNScoreShape(i, j2) << "\t"
				<< _nullWedge->GetNScorePN(i, j2) << "\n";
			}
		}
		
		ofs_nshrt.close();
	
		ofstream ofs_null(nullOut.c_str(), ofstream::out);
		
		ofs_null << "i\tj\tk\tStruct\tFreq\tDepthF\tShape\tN\n";
		
		for(int i = 0; i < _depth; i++)
		{
			if(i == 0)
			{
				ofs_null << i << "\t" << 0 << "\t" << 0 << "\t" 
				<< _nullWedge->GetScoreM(i, 0, 0) << "\t" 
				<< _nullWedge->GetFreqM(i, 0, 0) << "\t"
				<< _nullWedge->GetDepthFreqM(i, 0, 0) << "\t"
				<< _nullWedge->GetScoreD(i, 0, 0) << "\t" 
				<< _nullWedge->GetShapeN(i, 0, 0) << "\n";
			}
			
			for(int j = 0; j < i; j++)
			{
				if(j == 0)
				{
					ofs_null << i << "\t" << j << "\t" << 0 << "\t" 
					<< _nullWedge->GetScoreM(i, j, 0) << "\t" 
					<< _nullWedge->GetFreqM(i, j, 0) << "\t"
					<< _nullWedge->GetDepthFreqM(i, j, 0) << "\t"
					<< _nullWedge->GetScoreD(i, j, 0) << "\t" 
					<< _nullWedge->GetShapeN(i, j, 0) << "\n";
				}
				
				for(int l = 0; l < j; l++)
				{
					ofs_null << i << "\t" << j << "\t" << l << "\t" 
					<< _nullWedge->GetScoreM(i, j, l) << "\t" 
					<< _nullWedge->GetFreqM(i, j, l) << "\t"
					<< _nullWedge->GetDepthFreqM(i, j, l) << "\t"
					<< _nullWedge->GetScoreD(i, j, l) << "\t" 
					<< _nullWedge->GetShapeN(i, j, l) << "\n";
				}
			}
		}
		
		ofs_null.close();
	}
}

void TreeManager::SubtractAndResolve()
{
	if(_doNull)
	{
		_mainWedge->SubtractNull(false);
		_mainTree->ResolveParetos();
		
		cout << "\nMain Tree Output - POST-SUBTRACT\n";
		_mainTree->StateZSummary(false);
		
		//_mainWedge->MergeIn();
		//_nullWedge->MergeIn();
	}
	
	delete _mainTree;
}

void TreeManager::RunZDSCS(char * genome, bool hasDepths, char * depthFile)
{
	unsigned short tempMax = Node::NMAX;
	Node::NMAX = 0;
	
	_mainTree = new Tree(_depth, false, _doSubtract);
	_nullTree = new Tree(_depth, true, _doSubtract);
	
	if(_doSubtract)
	{
		cout << "Reading into First Tree" << endl;
		if(hasDepths)
		{
			_mainTree->ReadDepths(depthFile);
		}
		_mainTree->ReadFasta(genome);
		cout << "Summarising Tree, NMAX: " << Node::NMAX << endl;
		_mainTree->FirstTreeSummary();
	}
	
	double mainEscape = _mainTree->GetSEscape();
	delete _mainTree;
	
	Node::NMAX = tempMax;
	
	if(_doNull)
	{
		cout << "\nReading into Null Tree" << endl;
		_nullTree->SetSkip(mainEscape);
		if(hasDepths)
		{
			_nullTree->ReadDepths(depthFile);
		}
		_nullTree->ReadFasta(genome);
		cout << "Expanding Tree, NMAX: " << Node::NMAX << "\n" << endl;
		_nullTree->ExpandTree();
		cout << "Null Tree Output\n";
		_nullWedge = _nullTree->GetZWedge();

		_nullTree->ResolveParetos();
		_nullTree->StateZSummary(true);
	}
	
	delete _nullTree;
	
	_mainTree = new Tree(_depth, false, _doSubtract);
	
	cout << "\nReading into Main Tree" << endl;
	if(hasDepths)
	{
		_mainTree->ReadDepths(depthFile);
	}
	_mainTree->ReadFasta(genome);
	cout << "Expanding Tree, NMAX: " << Node::NMAX << "\n" << endl;
	_mainTree->ExpandTree();
	_mainWedge = _mainTree->GetZWedge();
	//_mainWedge->MergeIn();
	
	cout << "Main Tree Output - PRE-SUBTRACT\n";
	_mainTree->ResolveParetos();
	_mainTree->StateZSummary(true);
	if(_doNull)
	{
		_mainTree->SetNullReference(_nullWedge);
	}
}
