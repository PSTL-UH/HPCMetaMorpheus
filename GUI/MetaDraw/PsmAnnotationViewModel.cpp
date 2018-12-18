#include "PsmAnnotationViewModel.h"
#include "../../EngineLayer/MetaDraw/MetaDrawPsm.h"

using namespace OxyPlot;
using namespace OxyPlot::Axes;
using namespace OxyPlot::Series;
using namespace OxyPlot::Annotations;
using namespace MassSpectrometry;
using namespace Proteomics::Fragmentation;
using namespace EngineLayer;
using namespace Chemistry;

namespace ViewModels
{

std::unordered_map<ProductType*, OxyColor*> PsmAnnotationViewModel::productTypeDrawColors =
{
	{ProductType::b, OxyColors::Blue},
	{ProductType::y, OxyColors::Purple},
	{ProductType::c, OxyColors::Gold},
	{ProductType::zDot, OxyColors::Orange},
	{ProductType::D, OxyColors::DodgerBlue},
	{ProductType::M, OxyColors::Firebrick}
};
std::unordered_map<ProductType*, OxyColor*> PsmAnnotationViewModel::betaPeptideProductTypeDrawColors =
{
	{ProductType::b, OxyColors::LightBlue},
	{ProductType::y, OxyColors::MediumPurple},
	{ProductType::c, OxyColors::LightGoldenrodYellow},
	{ProductType::zDot, OxyColors::OrangeRed},
	{ProductType::D, OxyColors::AliceBlue},
	{ProductType::M, OxyColors::LightCoral}
};

	PlotModel *PsmAnnotationViewModel::getModel() const
	{
		return this->privateModel;
	}

	void PsmAnnotationViewModel::setModel(PlotModel *value)
	{
		this->privateModel = value;
		NotifyPropertyChanged(L"Model");
	}

	void PsmAnnotationViewModel::NotifyPropertyChanged(const std::wstring &propertyName)
	{
		PropertyChangedEventHandler handler = PropertyChanged;
		if (handler != nullptr)
		{
			PropertyChangedEventArgs tempVar(propertyName);
			handler(this, &tempVar);
		}
	}

	PsmAnnotationViewModel::PsmAnnotationViewModel()
	{
		// Create and Set the Model property, the INotifyPropertyChanged event will make the WPF Plot control update its content
		privateModel = new PlotModel();
		privateModel->Title = L"Spectrum Annotation";
		privateModel->Subtitle = L"using OxyPlot";
	}

	void PsmAnnotationViewModel::DrawPeptideSpectralMatch(MsDataScan *msDataScan, MetaDrawPsm *psmToDraw)
	{
		// x is m/z, y is intensity
		auto spectrumMzs = msDataScan->MassSpectrum.XArray;
		auto spectrumIntensities = msDataScan->MassSpectrum.YArray;

		std::wstring subtitle = psmToDraw->getFullSequence();
		if (psmToDraw != nullptr && psmToDraw->getBetaPeptideBaseSequence() != L"")
		{
			subtitle = psmToDraw->getFullSequence() + L"\n" + psmToDraw->getBetaPeptideFullSequence();
		}
		PlotModel *model = new PlotModel();
		model->Title = L"Spectrum Annotation of Scan #" + msDataScan->OneBasedScanNumber;
		model->DefaultFontSize = 15;
		model->Subtitle = subtitle;
		LinearAxis *tempVar = new LinearAxis();
		tempVar->Position = AxisPosition::Bottom;
		tempVar->Title = L"m/z";
		tempVar->Minimum = 0;
		tempVar->Maximum = spectrumMzs->Max() * 1.02;
		tempVar->AbsoluteMinimum = 0;
		tempVar->AbsoluteMaximum = spectrumMzs->Max() * 5;
		model->Axes->Add(tempVar);
		LinearAxis *tempVar2 = new LinearAxis();
		tempVar2->Position = AxisPosition->Left;
		tempVar2->Title = L"Intensity";
		tempVar2->Minimum = 0;
		tempVar2->Maximum = spectrumIntensities->Max() * 1.2;
		tempVar2->AbsoluteMinimum = 0;
		tempVar2->AbsoluteMaximum = spectrumIntensities->Max() * 1.3;
		model->Axes->Add(tempVar2);
		model->Axes[1].Zoom(0, spectrumIntensities->Max() * 1.1);

		std::vector<LineSeries*> allIons(spectrumMzs->Length);

		// draw the matched peaks; if the PSM is null, we're just drawing the peaks in the scan without annotation, so skip this part
		if (psmToDraw != nullptr)
		{
			for (auto peak : psmToDraw->getMatchedIons())
			{
				OxyColor *ionColor;

				if (productTypeDrawColors.find(peak->NeutralTheoreticalProduct.ProductType) != productTypeDrawColors.end())
				{
					ionColor = productTypeDrawColors[peak->NeutralTheoreticalProduct.ProductType];
				}
				else
				{
					ionColor = OxyColors::Turquoise;
				}

				int i = msDataScan->MassSpectrum.GetClosestPeakIndex(peak->NeutralTheoreticalProduct.NeutralMass.ToMz(peak->Charge))->Value;

				// peak line
				allIons[i] = new LineSeries();
				allIons[i]->Color = ionColor;
				allIons[i]->StrokeThickness = STROKE_THICKNESS_ANNOTATED;
				DataPoint tempVar3(peak->Mz, 0);
				allIons[i]->Points->Add(&tempVar3);
				DataPoint tempVar4(peak->Mz, spectrumIntensities[i]);
				allIons[i]->Points->Add(&tempVar4);

				// peak annotation
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
				std::wstring peakAnnotationText = peak->NeutralTheoreticalProduct.ProductType.ToString()->ToLower() + peak->NeutralTheoreticalProduct.TerminusFragment.FragmentNumber + L" (" + peak->Mz.ToString(L"F3") + L")";
				if (peak->NeutralTheoreticalProduct.NeutralLoss != 0)
				{
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
					peakAnnotationText = peak->NeutralTheoreticalProduct.ProductType.ToString()->ToLower() + peak->NeutralTheoreticalProduct.TerminusFragment.FragmentNumber + L"-" + peak->NeutralTheoreticalProduct.NeutralLoss.ToString(L"F2") + L" (" + peak->Mz.ToString(L"F3") + L")";
				}

				auto peakAnnotation = new TextAnnotation();
				peakAnnotation->TextRotation = -60;
				peakAnnotation->Font = L"Arial";
				peakAnnotation->FontSize = 12;
				peakAnnotation->FontWeight = 2.0;
				peakAnnotation->TextColor = ionColor;
				peakAnnotation->StrokeThickness = 0;
				peakAnnotation->Text = peakAnnotationText;
				peakAnnotation->TextPosition = new DataPoint(allIons[i]->Points[1].X, allIons[i]->Points[1].Y + peakAnnotation->Text->Length * 1.5 / 4);
				peakAnnotation->TextHorizontalAlignment = HorizontalAlignment->Left;
				model->Annotations->Add(peakAnnotation);

				model->Series->Add(allIons[i]);

//C# TO C++ CONVERTER TODO TASK: A 'delete peakAnnotation' statement was not added since peakAnnotation was passed to a method or constructor. Handle memory management manually.
			}

			if (psmToDraw->getBetaPeptideBaseSequence() != L"")
			{
				for (auto peak : psmToDraw->getBetaPeptideMatchedIons())
				{
					OxyColor *ionColor;

					if (productTypeDrawColors.find(peak->NeutralTheoreticalProduct.ProductType) != productTypeDrawColors.end())
					{
						ionColor = betaPeptideProductTypeDrawColors[peak->NeutralTheoreticalProduct.ProductType];
					}
					else
					{
						ionColor = OxyColors::Turquoise;
					}

					int i = msDataScan->MassSpectrum.GetClosestPeakIndex(peak->NeutralTheoreticalProduct.NeutralMass.ToMz(peak->Charge))->Value;

					// peak line
					allIons[i] = new LineSeries();
					allIons[i]->Color = ionColor;
					allIons[i]->StrokeThickness = STROKE_THICKNESS_ANNOTATED;
					DataPoint tempVar5(peak->Mz, 0);
					allIons[i]->Points->Add(&tempVar5);
					DataPoint tempVar6(peak->Mz, spectrumIntensities[i]);
					allIons[i]->Points->Add(&tempVar6);

					// peak annotation
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
					std::wstring peakAnnotationText = L"beta-" + peak->NeutralTheoreticalProduct.ProductType.ToString()->ToLower() + peak->NeutralTheoreticalProduct.TerminusFragment.FragmentNumber + L" (" + peak->Mz.ToString(L"F3") + L")";
					if (peak->NeutralTheoreticalProduct.NeutralLoss != 0)
					{
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
						peakAnnotationText = L"beta-" + peak->NeutralTheoreticalProduct.ProductType.ToString()->ToLower() + peak->NeutralTheoreticalProduct.TerminusFragment.FragmentNumber + L"-" + peak->NeutralTheoreticalProduct.NeutralLoss.ToString(L"F2") + L" (" + peak->Mz.ToString(L"F3") + L")";
					}

					auto peakAnnotation = new TextAnnotation();
					peakAnnotation->TextRotation = -60;
					peakAnnotation->Font = L"Arial";
					peakAnnotation->FontSize = 12;
					peakAnnotation->FontWeight = 2.0;
					peakAnnotation->TextColor = ionColor;
					peakAnnotation->StrokeThickness = 0;
					peakAnnotation->Text = peakAnnotationText;
					peakAnnotation->TextPosition = new DataPoint(allIons[i]->Points[1].X, allIons[i]->Points[1].Y + peakAnnotation->Text->Length * 1.5 / 4);
					peakAnnotation->TextHorizontalAlignment = HorizontalAlignment->Left;
					model->Annotations->Add(peakAnnotation);

					model->Series->Add(allIons[i]);

//C# TO C++ CONVERTER TODO TASK: A 'delete peakAnnotation' statement was not added since peakAnnotation was passed to a method or constructor. Handle memory management manually.
				}
			}
		}

		// draw the remaining unmatched peaks
		for (int i = 0; i < spectrumMzs->Length; i++)
		{
			// peak has already been drawn (it is a matched peak)
			if (allIons[i] != nullptr)
			{
				continue;
			}

			allIons[i] = new LineSeries();
			allIons[i]->Color = OxyColors::DimGray;
			allIons[i]->StrokeThickness = STROKE_THICKNESS_UNANNOTATED;
			DataPoint tempVar7(spectrumMzs[i], 0);
			allIons[i]->Points->Add(&tempVar7);
			DataPoint tempVar8(spectrumMzs[i], spectrumIntensities[i]);
			allIons[i]->Points->Add(&tempVar8);
			model->Series->Add(allIons[i]);
		}

		// Axes are created automatically if they are not defined

		// Set the Model property, the INotifyPropertyChanged event will make the WPF Plot control update its content
		this->setModel(model);

//C# TO C++ CONVERTER TODO TASK: A 'delete model' statement was not added since model was assigned to another object. Handle memory management manually.
	}
}
