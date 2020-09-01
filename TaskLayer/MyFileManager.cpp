#include "MyFileManager.h"
#include "../EngineLayer/CommonParameters.h"
#include "../EngineLayer/EventArgs/StringEventArgs.h"
#include "MetaMorpheusTask.h"

#include "MzML/Mzml.h"
using namespace IO::MzML;

#include "Mgf/Mgf.h"
using namespace IO::Mgf;

using namespace EngineLayer;

#include "MassSpectrometry/MassSpectrometry.h"
using namespace MassSpectrometry;


namespace TaskLayer
{

    MyFileManager::MyFileManager(bool disposeOfFileWhenDone) : DisposeOfFileWhenDone(disposeOfFileWhenDone)
    {
    }
    
    
    MsDataFile *MyFileManager::LoadFile(const std::string &origDataFile, std::optional<int> topNpeaks,
                                        std::optional<double> minRatio, bool trimMs1Peaks, bool trimMsMsPeaks,
                                        CommonParameters *commonParameters)
    {

        std::unordered_map<std::string, MsDataFile*>::const_iterator MyMsDataFiles_iterator = MyMsDataFiles.find(origDataFile);
        if (MyMsDataFiles_iterator != MyMsDataFiles.end() && MyMsDataFiles_iterator->second != nullptr )
        {
            return MyMsDataFiles_iterator->second;            
        }
        
        // By now know that need to load this file!!!
        //std::lock_guard<std::mutex> lock(FileLoadingLock);
        std::string extension = origDataFile.substr(origDataFile.find_last_of("."));
        
        if ( extension == ".mzML" )
        {
            FilteringParams *filter = new FilteringParams(topNpeaks, minRatio, 1, trimMs1Peaks, trimMsMsPeaks);
            MyMsDataFiles[origDataFile] = Mzml::LoadAllStaticData(origDataFile, filter,
                                                                  commonParameters->getMaxThreadsToUsePerFile());
            delete filter;        
        }
        else if ( extension == ".mgf" )
        {
            // mgf reader does not handle filter as of now.
            MyMsDataFiles[origDataFile] = Mgf::LoadAllStaticData(origDataFile);
        }
        else
        {
            Warn("No capability for reading " + origDataFile);
        }
        
        return MyMsDataFiles[origDataFile];
    }

    void MyFileManager::DoneWithFile(const std::string &origDataFile)
    {
        if (DisposeOfFileWhenDone)
        {
            MyMsDataFiles[origDataFile] = nullptr;
        }
    }
    
    void MyFileManager::Warn(const std::string &v)
    {
        MetaMorpheusTask::Warn(v);
    }
}
