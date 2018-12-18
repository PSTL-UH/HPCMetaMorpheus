#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <tuple>

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class PeptideSpectralMatch; }

using namespace Proteomics;
using namespace Proteomics::AminoAcidPolymer;

namespace EngineLayer
{
	namespace HistogramAnalysis
	{
		class Bin
		{
		private:
			int privatePepNlocCount = 0;
			int privatePepClocCount = 0;
			int privateProtNlocCount = 0;
			int privateProtClocCount = 0;
			std::wstring privateCombos = L"-";
			std::wstring privateUnimodDiffs = L"-";
			std::wstring privateUniprotID = L"-";
			std::wstring privateUnimodFormulas = L"-";
			std::wstring privateUnimodId = L"-";
			double privateMassShift = 0;
			std::wstring privateMine;
			std::unordered_map<wchar_t, int> privateAAsInCommon;
			int privateOverlapping = 0;
			double privateFracWithSingle = 0;
			double privateMedianLength = 0;

		public:
			std::wstring AA = L"-";
			std::unordered_map<wchar_t, int> ResidueCount;
			std::unordered_map<std::wstring, std::tuple<std::wstring, std::wstring, PeptideSpectralMatch*>> UniquePSMs;
			std::unordered_map<std::wstring, int> ModsInCommon;

			Bin(double massShift);

				int getPepNlocCount() const;
				void setPepNlocCount(int value);
				int getPepClocCount() const;
				void setPepClocCount(int value);
				int getProtNlocCount() const;
				void setProtNlocCount(int value);
				int getProtClocCount() const;
				void setProtClocCount(int value);
				std::wstring getCombos() const;
				void setCombos(const std::wstring &value);
				std::wstring getUnimodDiffs() const;
				void setUnimodDiffs(const std::wstring &value);
				std::wstring getUniprotID() const;
				void setUniprotID(const std::wstring &value);
				std::wstring getUnimodFormulas() const;
				void setUnimodFormulas(const std::wstring &value);
				std::wstring getUnimodId() const;
				void setUnimodId(const std::wstring &value);
				double getMassShift() const;

			int getCount() const;

			int getCountDecoy() const;

			int getCountTarget() const;

			int getLocalizeableTarget() const;

				std::wstring getMine() const;
				void setMine(const std::wstring &value);
				std::unordered_map<wchar_t, int> getAAsInCommon() const;
				void setAAsInCommon(const std::unordered_map<wchar_t, int> &value);
				int getOverlapping() const;
				void setOverlapping(int value);
				double getFracWithSingle() const;
				void setFracWithSingle(double value);
				double getMedianLength() const;
				void setMedianLength(double value);

			void IdentifyResidues();

			void IdentifyCombos(double v, std::unordered_set<std::tuple<double, double, double>> &ok);

			double ComputeZ(double v);

			void IdentifyUniprotBins(double v);

			void IdentifyAA(double v);

			void IdentifyUnimodBins(double v);

			void Add(PeptideSpectralMatch *ok);
		};
	}
}
