#pragma once

namespace TaskLayer
{
	class CalibrationParameters
	{
	private:
		bool privateWriteIntermediateFiles = false;
		int privateMinMS1IsotopicPeaksNeededForConfirmedIdentification = 0;
		int privateMinMS2IsotopicPeaksNeededForConfirmedIdentification = 0;
		int privateNumFragmentsNeededForEveryIdentification = 0;

	public:
		CalibrationParameters();

		bool getWriteIntermediateFiles() const;
		void setWriteIntermediateFiles(bool value);

		int getMinMS1IsotopicPeaksNeededForConfirmedIdentification() const;
		void setMinMS1IsotopicPeaksNeededForConfirmedIdentification(int value);
		int getMinMS2IsotopicPeaksNeededForConfirmedIdentification() const;
		void setMinMS2IsotopicPeaksNeededForConfirmedIdentification(int value);
		int getNumFragmentsNeededForEveryIdentification() const;
		void setNumFragmentsNeededForEveryIdentification(int value);
	};
}
