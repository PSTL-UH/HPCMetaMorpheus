#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include "tangible_event.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class MetaDrawPsm; }

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
	class PsmAnnotationViewModel : public INotifyPropertyChanged
	{
	private:
		static constexpr double STROKE_THICKNESS_UNANNOTATED = 0.5;
		static constexpr double STROKE_THICKNESS_ANNOTATED = 2.0;
		PlotModel *privateModel;

		static std::unordered_map<ProductType*, OxyColor*> productTypeDrawColors;

		static std::unordered_map<ProductType*, OxyColor*> betaPeptideProductTypeDrawColors;

	public:
		virtual ~PsmAnnotationViewModel()
		{
			delete privateModel;
		}

		PlotModel *getModel() const;
		void setModel(PlotModel *value);

		TangibleEvent<PropertyChangedEventHandler> *PropertyChanged = new TangibleEvent<PropertyChangedEventHandler>();

	protected:
		void NotifyPropertyChanged(const std::wstring &propertyName);

	public:
		PsmAnnotationViewModel();

		// single peptides (not crosslink)
		void DrawPeptideSpectralMatch(MsDataScan *msDataScan, MetaDrawPsm *psmToDraw);

	};
}
