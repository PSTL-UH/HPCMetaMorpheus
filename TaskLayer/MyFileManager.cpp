#include "MyFileManager.h"
#include "../EngineLayer/CommonParameters.h"
#include "../EngineLayer/EventArgs/StringEventArgs.h"
#include <experimental/filesystem>

#include "MzML/Mzml.h"
using namespace IO::MzML;

using namespace EngineLayer;
//using namespace IO::Mgf;
#if defined(NETFRAMEWORK)
using namespace IO::Thermo;
#endif

#include "MassSpectrometry/MassSpectrometry.h"
using namespace MassSpectrometry;

namespace TaskLayer
{

const std::string MyFileManager::AssumedThermoMsFileReaderDllPath = "Thermo/MSFileReader";
const std::string MyFileManager::DesiredFileIoVersion = "3.0";
const std::string MyFileManager::DesiredFregistryVersion = "3.0";
const std::string MyFileManager::DesiredXRawFileVersion = "3.0.29.0";

    MyFileManager::MyFileManager(bool disposeOfFileWhenDone) : DisposeOfFileWhenDone(disposeOfFileWhenDone)
    {
        WarnHandler = new EventHandler<StringEventArgs>();
    }
    
    bool MyFileManager::SeeIfOpen(std::string path)
    {
        return ( MyMsDataFiles.find(path) != MyMsDataFiles.end() && 
                 MyMsDataFiles[path] != nullptr ); 
    }
    
    MyFileManager::ThermoMsFileReaderVersionCheck MyFileManager::ValidateThermoMsFileReaderVersion()
    {
        std::string fileIoAssumedPath = AssumedThermoMsFileReaderDllPath + "Fileio_x64.dll";
        std::string fregistryAssumedPath = AssumedThermoMsFileReaderDllPath + "fregistry_x64.dll";
        std::string xRawFileAssumedPath = AssumedThermoMsFileReaderDllPath + "XRawfile2_x64.dll";
        
        if ( std::experimental::filesystem::exists(fileIoAssumedPath)    &&
             std::experimental::filesystem::exists(fregistryAssumedPath) &&
             std::experimental::filesystem::exists(xRawFileAssumedPath))
        {
            //std::string fileIoVersion = FileVersionInfo::GetVersionInfo(fileIoAssumedPath)->FileVersion;
            //std::string fregistryVersion = FileVersionInfo::GetVersionInfo(fregistryAssumedPath)->FileVersion;
            //std::string xRawFileVersion = FileVersionInfo::GetVersionInfo(xRawFileAssumedPath)->FileVersion;
            
            // EDGAR: for now just set the file version to what is required.
            std::string fileIoVersion = "3.0";
            std::string fregistryVersion = "3.0";
            std::string xRawFileVersion = "3.0.29.0";
            
            if ( fileIoVersion == DesiredFileIoVersion       &&
                 fregistryVersion == DesiredFregistryVersion &&
                 xRawFileVersion == DesiredXRawFileVersion)
            {
                return ThermoMsFileReaderVersionCheck::CorrectVersion;
            }
            else
            {
                return ThermoMsFileReaderVersionCheck::IncorrectVersion;
            }
        }
        else if ( std::experimental::filesystem::exists(fileIoAssumedPath)    ||
                  std::experimental::filesystem::exists(fregistryAssumedPath) ||
                  std::experimental::filesystem::exists(xRawFileAssumedPath))
        {
            return ThermoMsFileReaderVersionCheck::SomeDllsMissing;
        }
        
        return ThermoMsFileReaderVersionCheck::DllsNotFound;
    }
    
    MsDataFile *MyFileManager::LoadFile(const std::string &origDataFile, std::optional<int> topNpeaks,
                                        std::optional<double> minRatio, bool trimMs1Peaks, bool trimMsMsPeaks,
                                        CommonParameters *commonParameters)
    {
        FilteringParams *filter = new FilteringParams(topNpeaks, minRatio, 1, trimMs1Peaks, trimMsMsPeaks);
        //MsDataFile value = new MsDataFile(false);
        std::unordered_map<std::string, MsDataFile*>::const_iterator MyMsDataFiles_iterator = MyMsDataFiles.find(origDataFile);
        if (MyMsDataFiles_iterator != MyMsDataFiles.end() )
        {
            delete filter;
            return MyMsDataFiles_iterator->second;
            
        }
        
        // By now know that need to load this file!!!
        //std::lock_guard<std::mutex> lock(FileLoadingLock);
        std::string extension = origDataFile.substr(origDataFile.find_last_of("."));
        
        if ( extension == ".mzM" )
        {
            // arguments:  const std::string &filePath, FilteringParams *filterParams, int maxThreads
            MyMsDataFiles[origDataFile] = Mzml::LoadAllStaticData(origDataFile, filter,
                                                                  commonParameters->getMaxThreadsToUsePerFile());
        }
#ifdef ORIG
        // MGF files are not support right now.
        else if ( extension == ".mgf" )
        {
            MyMsDataFiles[origDataFile] = Mgf::LoadAllStaticData(origDataFile, filter);
        }
#endif
        else
        {
#if defined(NETFRAMEWORK)
            MyMsDataFiles[origDataFile] = ThermoStaticData::LoadAllStaticData(origDataFile, filter);
#else
            Warn("No capability for reading " + origDataFile);
#endif
        }
        
        //C# TO C++ CONVERTER TODO TASK: A 'delete filter' statement was not added since filter was 
        //passed to a method or constructor. Handle memory management manually.
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
        std::vector<std::string> tmpvec;
        StringEventArgs tempVar(v, tmpvec);
        if ( WarnHandler != nullptr ) {
            WarnHandler->Invoke(tempVar);
        }
    }
}
