#pragma once

#include <string>
#include <unordered_map>
#include <optional>
#include <mutex>
#include "tangible_filesystem.h"
#include "EventHandler.h"
#include "../EngineLayer/EventArgs/StringEventArgs.h"

#include "../EngineLayer/CommonParameters.h"
using namespace EngineLayer;


using namespace MassSpectrometry;

namespace TaskLayer
{
    class MyFileManager
    {
    public:
        enum class ThermoMsFileReaderVersionCheck
        {
            DllsNotFound,
            IncorrectVersion,
            CorrectVersion,
            SomeDllsMissing
        };
        
    private:
        const bool DisposeOfFileWhenDone;
        std::unordered_map<std::string, MsDataFile*> MyMsDataFiles = std::unordered_map<std::string, MsDataFile*>();
        std::mutex FileLoadingLock;
        
    public:
        MyFileManager(bool disposeOfFileWhenDone);
        
       
        bool SeeIfOpen(const std::string path);
                
        MsDataFile *LoadFile(const std::string &origDataFile, std::optional<int> topNpeaks, std::optional<double> minRatio,
                             bool trimMs1Peaks, bool trimMsMsPeaks, CommonParameters *commonParameters);
        
        void DoneWithFile(const std::string &origDataFile);
        
    private:
        void Warn(const std::string &v);
    };
}
