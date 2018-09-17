##
## Auto Generated makefile by CodeLite IDE
## any manual changes will be erased      
##
## Release
ProjectName            :=UGPep
ConfigurationName      :=Release
WorkspacePath          :=.
ProjectPath            :=.
IntermediateDirectory  :=./Release
OutDir                 := $(IntermediateDirectory)
CurrentFileName        :=
CurrentFilePath        :=
CurrentFileFullPath    :=
User                   :=Oliver Rimington
Date                   :=07/05/18
CodeLitePath           :=
LinkerName             :=/usr/bin/g++
SharedObjectLinkerName :=/usr/bin/g++ -shared -fPIC
ObjectSuffix           :=.o
DependSuffix           :=.o.d
PreprocessSuffix       :=.i
DebugSwitch            :=-g 
IncludeSwitch          :=-I
LibrarySwitch          :=-l
OutputSwitch           :=-o 
LibraryPathSwitch      :=-L
PreprocessorSwitch     :=-D
SourceSwitch           :=-c 
OutputFile             :=$(IntermediateDirectory)/$(ProjectName)
Preprocessors          :=$(PreprocessorSwitch)NDEBUG 
ObjectSwitch           :=-o 
ArchiveOutputSwitch    := 
PreprocessOnlySwitch   :=-E
ObjectsFileList        :="UGPep.txt"
PCHCompileFlags        :=
MakeDirCommand         :=mkdir -p
LinkOptions            :=  -pthread
IncludePath            :=  $(IncludeSwitch). $(IncludeSwitch). 
IncludePCH             := 
RcIncludePath          := 
Libs                   := 
ArLibs                 :=  
LibPath                := $(LibraryPathSwitch). 

##
## Common variables
## AR, CXX, CC, AS, CXXFLAGS and CFLAGS can be overriden using an environment variables
##
AR       := /usr/bin/ar rcu
CXX      := /usr/bin/g++
CC       := /usr/bin/gcc
CXXFLAGS :=  -O2 -Wall -std=c++11  $(Preprocessors)
CFLAGS   :=  -O2 -Wall $(Preprocessors)
ASFLAGS  := 
AS       := /usr/bin/as


##
## User defined environment variables
##
Objects0=$(IntermediateDirectory)/main.cpp$(ObjectSuffix) $(IntermediateDirectory)/Node.cpp$(ObjectSuffix) $(IntermediateDirectory)/Tree.cpp$(ObjectSuffix) $(IntermediateDirectory)/WedgeMatrix.cpp$(ObjectSuffix) $(IntermediateDirectory)/TreeManager.cpp$(ObjectSuffix) $(IntermediateDirectory)/ThreadManager.cpp$(ObjectSuffix) 



Objects=$(Objects0) 

##
## Main Build Targets 
##
.PHONY: all clean PreBuild PrePreBuild PostBuild MakeIntermediateDirs
all: $(OutputFile)

$(OutputFile): $(IntermediateDirectory)/.d $(Objects) 
	@$(MakeDirCommand) $(@D)
	@echo "" > $(IntermediateDirectory)/.d
	@echo $(Objects0)  > $(ObjectsFileList)
	$(LinkerName) $(OutputSwitch)$(OutputFile) @$(ObjectsFileList) $(LibPath) $(Libs) $(LinkOptions)

MakeIntermediateDirs:
	@test -d ./Release || $(MakeDirCommand) ./Release


$(IntermediateDirectory)/.d:
	@test -d ./Release || $(MakeDirCommand) ./Release

PreBuild:


##
## Objects
##
$(IntermediateDirectory)/main.cpp$(ObjectSuffix): main.cpp $(IntermediateDirectory)/main.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "main.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/main.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/main.cpp$(DependSuffix): main.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/main.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/main.cpp$(DependSuffix) -MM "main.cpp"

$(IntermediateDirectory)/main.cpp$(PreprocessSuffix): main.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/main.cpp$(PreprocessSuffix) "main.cpp"

$(IntermediateDirectory)/Node.cpp$(ObjectSuffix): Node.cpp $(IntermediateDirectory)/Node.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "Node.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/Node.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/Node.cpp$(DependSuffix): Node.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/Node.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/Node.cpp$(DependSuffix) -MM "Node.cpp"

$(IntermediateDirectory)/Node.cpp$(PreprocessSuffix): Node.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/Node.cpp$(PreprocessSuffix) "Node.cpp"

$(IntermediateDirectory)/Tree.cpp$(ObjectSuffix): Tree.cpp $(IntermediateDirectory)/Tree.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "Tree.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/Tree.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/Tree.cpp$(DependSuffix): Tree.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/Tree.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/Tree.cpp$(DependSuffix) -MM "Tree.cpp"

$(IntermediateDirectory)/Tree.cpp$(PreprocessSuffix): Tree.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/Tree.cpp$(PreprocessSuffix) "Tree.cpp"


$(IntermediateDirectory)/WedgeMatrix.cpp$(ObjectSuffix): WedgeMatrix.cpp $(IntermediateDirectory)/WedgeMatrix.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "WedgeMatrix.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/WedgeMatrix.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/WedgeMatrix.cpp$(DependSuffix): WedgeMatrix.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/WedgeMatrix.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/WedgeMatrix.cpp$(DependSuffix) -MM "WedgeMatrix.cpp"

$(IntermediateDirectory)/WedgeMatrix.cpp$(PreprocessSuffix): WedgeMatrix.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/WedgeMatrix.cpp$(PreprocessSuffix) "WedgeMatrix.cpp"

$(IntermediateDirectory)/TreeManager.cpp$(ObjectSuffix): TreeManager.cpp $(IntermediateDirectory)/TreeManager.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "TreeManager.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/TreeManager.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/TreeManager.cpp$(DependSuffix): TreeManager.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/TreeManager.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/TreeManager.cpp$(DependSuffix) -MM "TreeManager.cpp"

$(IntermediateDirectory)/TreeManager.cpp$(PreprocessSuffix): TreeManager.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/TreeManager.cpp$(PreprocessSuffix) "TreeManager.cpp"

$(IntermediateDirectory)/ThreadManager.cpp$(ObjectSuffix): ThreadManager.cpp $(IntermediateDirectory)/ThreadManager.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "ThreadManager.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/ThreadManager.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/ThreadManager.cpp$(DependSuffix): ThreadManager.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/ThreadManager.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/ThreadManager.cpp$(DependSuffix) -MM "ThreadManager.cpp"

$(IntermediateDirectory)/ThreadManager.cpp$(PreprocessSuffix): ThreadManager.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/ThreadManager.cpp$(PreprocessSuffix) "ThreadManager.cpp"


-include $(IntermediateDirectory)/*$(DependSuffix)
##
## Clean
##
clean:
	$(RM) -r ./Release/


