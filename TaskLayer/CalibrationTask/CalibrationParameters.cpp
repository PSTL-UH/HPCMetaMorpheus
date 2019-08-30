#include "CalibrationParameters.h"

namespace TaskLayer
{

	CalibrationParameters::CalibrationParameters()
	{
		setWriteIntermediateFiles(false);
		setMinMS1IsotopicPeaksNeededForConfirmedIdentification(3);
		setMinMS2IsotopicPeaksNeededForConfirmedIdentification(2);
		setNumFragmentsNeededForEveryIdentification(10);
	}

	bool CalibrationParameters::getWriteIntermediateFiles() const
	{
		return privateWriteIntermediateFiles;
	}

	void CalibrationParameters::setWriteIntermediateFiles(bool value)
	{
		privateWriteIntermediateFiles = value;
	}

	int CalibrationParameters::getMinMS1IsotopicPeaksNeededForConfirmedIdentification() const
	{
		return privateMinMS1IsotopicPeaksNeededForConfirmedIdentification;
	}

	void CalibrationParameters::setMinMS1IsotopicPeaksNeededForConfirmedIdentification(int value)
	{
		privateMinMS1IsotopicPeaksNeededForConfirmedIdentification = value;
	}

	int CalibrationParameters::getMinMS2IsotopicPeaksNeededForConfirmedIdentification() const
	{
		return privateMinMS2IsotopicPeaksNeededForConfirmedIdentification;
	}

	void CalibrationParameters::setMinMS2IsotopicPeaksNeededForConfirmedIdentification(int value)
	{
		privateMinMS2IsotopicPeaksNeededForConfirmedIdentification = value;
	}

	int CalibrationParameters::getNumFragmentsNeededForEveryIdentification() const
	{
		return privateNumFragmentsNeededForEveryIdentification;
	}

	void CalibrationParameters::setNumFragmentsNeededForEveryIdentification(int value)
	{
		privateNumFragmentsNeededForEveryIdentification = value;
	}
}
