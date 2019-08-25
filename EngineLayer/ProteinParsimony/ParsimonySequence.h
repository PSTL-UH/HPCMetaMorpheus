#pragma once

#include <string>
#include <any>

using namespace Proteomics::ProteolyticDigestion;

namespace EngineLayer
{
	namespace ProteinParsimony
	{
		class ParsimonySequence
		{
		private:
			std::string privateSequence;
			Protease *privateProtease;

		public:
			ParsimonySequence(PeptideWithSetModifications *pwsm, bool TreatModPeptidesAsDifferentPeptides);

                        std::string getSequence() const;
                        Protease *getProtease() const;

			bool Equals(std::any obj);

			int GetHashCode();
		};
	}
}
