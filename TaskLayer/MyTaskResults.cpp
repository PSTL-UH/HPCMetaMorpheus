#include "MyTaskResults.h"
#include "DbForTask.h"
//#include "MetaMorpheusTask.h"
#include "stringbuilder.h"
#include "stringhelper.h"

namespace TaskLayer
{

    MyTaskResults::MyTaskResults(MetaMorpheusTask *s) : resultTexts(std::vector<std::string>())
    {
    }
    
    std::string MyTaskResults::ToString()
    {
        StringBuilder *sb = new StringBuilder();
        sb->appendLine("Time to run task: " + std::to_string(Time));
        sb->appendLine();
        sb->appendLine();
        sb->appendLine("--------------------------------------------------");
        //if ((NewSpectra.size() > 0 && NewSpectra.Any()) || (NewDatabases.size() > 0 && NewDatabases.Any()))
        if ( NewSpectra.size() > 0 || NewDatabases.size() )
        {
            sb->appendLine();
            sb->appendLine();
            sb->appendLine("New files:");
            if (NewSpectra.size() > 0 )
            {
                sb->appendLine("New spectra: ");
                sb->appendLine();
                //sb->appendLine(std::string::Join("\r\n" + "\t", NewSpectra));
                std::string delimiter = "\r\n\t";
                sb->appendLine(StringHelper::join(NewSpectra, delimiter));
            }
            if (NewDatabases.size() > 0 )
            {
                sb->appendLine("New databases: ");
#ifdef ORIG
                sb->appendLine(std::string::Join("\r\n" + "\t", NewDatabases.Select([&] (std::any b)   {
                                b::FilePath;
                            })).ToString());
#endif
                std::vector<std::string> vs;
                for ( auto b: NewDatabases ) {
                    vs.push_back(b->getFilePath());
                }
                std::string delimiter = "\r\n\t";
                sb->appendLine(StringHelper::join(vs, delimiter));
            }
            sb->appendLine();
            sb->appendLine();
            sb->appendLine("--------------------------------------------------");
        }
        sb->appendLine();
        sb->appendLine();
        sb->appendLine(niceText->toString());
        sb->appendLine();
        sb->appendLine();
        sb->appendLine("--------------------------------------------------");
        sb->appendLine();
        sb->appendLine();
        sb->appendLine("Engine Results:");
        sb->appendLine();
        for (auto ok : resultTexts)
        {
            sb->appendLine(ok);
            sb->appendLine();
        }
        
        std::string s =  sb->toString();
        delete sb;
        return s;
    }
    
    void MyTaskResults::AddResultText(const std::string &resultsText)
    {
        resultTexts.push_back(resultsText);
    }
    
    void MyTaskResults::AddNiceText(const std::string &niceTextString)
    {
        niceText->appendLine(niceTextString);
    }
}
