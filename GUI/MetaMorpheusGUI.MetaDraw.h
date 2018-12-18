#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include <stdexcept>
#include <any>
#include <functional>
#include <optional>
#include "stringhelper.h"
#include "tangible_filesystem.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace ViewModels { class PsmAnnotationViewModel; }
namespace TaskLayer { class MyFileManager; }
namespace EngineLayer { class MetaDrawPsm; }

using namespace ViewModels;
using namespace EngineLayer;
using namespace MassSpectrometry;
using namespace TaskLayer;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;

namespace MetaMorpheusGUI
{
	/// <summary>
	/// Interaction logic for MetaDraw.xaml
	/// </summary>
	class MetaDraw : public Window
	{
	private:
		PsmAnnotationViewModel *mainViewModel;
		MyFileManager *spectraFileManager;
		MsDataFile *MsDataFile;
		ObservableCollection<MetaDrawPsm*> *const peptideSpectralMatches;
		ICollectionView *peptideSpectralMatchesView;
		DataTable *const propertyView;
		std::wstring spectraFilePath;
		std::wstring tsvResultsFilePath;
		std::unordered_map<ProductType*, double> productTypeToYOffset;
		std::unordered_map<ProductType*, Color*> productTypeToColor;
		SolidColorBrush *modificationAnnotationColor;
		Regex *illegalInFileName = new Regex(LR"([\\/:*?"<>|])");

	public:
		virtual ~MetaDraw()
		{
			delete mainViewModel;
			delete spectraFileManager;
			delete MsDataFile;
			delete peptideSpectralMatches;
			delete peptideSpectralMatchesView;
			delete propertyView;
			delete modificationAnnotationColor;
			delete illegalInFileName;
		}

		MetaDraw();

	private:
		void SetUpDictionaries();

		void Window_Drop(std::any sender, DragEventArgs *e);

		void LoadFile(const std::wstring &filePath);

		void LoadPsms(const std::wstring &filename);

		void DrawPsm(int oneBasedScanNumber, const std::wstring &fullSequence = L"");

		/// <summary>
		/// Event triggers when a different cell is selected in the PSM data grid
		/// </summary>
		void dataGridScanNums_SelectedCellsChanged(std::any sender, SelectedCellsChangedEventArgs *e);

		void selectSpectraFileButton_Click(std::any sender, RoutedEventArgs *e);

		void selectPsmFileButton_Click(std::any sender, RoutedEventArgs *e);

//C# TO C++ CONVERTER TODO TASK: There is no equivalent in C++ to the 'async' keyword:
//ORIGINAL LINE: private async void loadFilesButton_Click(object sender, RoutedEventArgs e)
		void loadFilesButton_Click(std::any sender, RoutedEventArgs *e);

		void TextBox_TextChanged(std::any sender, TextChangedEventArgs *e);

		void DrawAnnotatedBaseSequence(MetaDrawPsm *psm);

		void dataGridProperties_SelectedCellsChanged(std::any sender, SelectedCellsChangedEventArgs *e);

		//private void PDFButton_Click(object sender, RoutedEventArgs e)
		//{
		//    if (dataGridScanNums.SelectedCells.Count == 0)
		//    {
		//        MessageBox.Show("Please select at least one scan to export");
		//    }

		//    int num = dataGridScanNums.SelectedItems.Count;
		//    string writeDirectory = Path.Combine(Directory.GetParent(tsvResultsFilePath).FullName, "PDF");

		//    if (!Directory.Exists(writeDirectory))
		//    {
		//        Directory.CreateDirectory(writeDirectory);
		//    }

		//    foreach (object selectedItem in dataGridScanNums.SelectedItems)
		//    {
		//        MetaDrawPsm psm = (MetaDrawPsm)selectedItem;
		//        string myString = illegalInFileName.Replace(psm.FullSequence, "").Substring(0, 40);
		//        ExportToPdf(psm, Path.Combine(writeDirectory, psm.Ms2ScanNumber + "_" + myString + ".pdf"));
		//    }

		//    dataGridScanNums.SelectedItem = dataGridScanNums.SelectedItem;
		//    MessageBox.Show(string.Format("{0} PDFs exported to " + writeDirectory, num));
		//}

		//private void ExportToPdf(MetaDrawPsm psm, string path)
		//{
		//    System.Reflection.PropertyInfo[] temp = psm.GetType().GetProperties();

		//    for (int i = 4; i < temp.Length; i++)
		//    {
		//        propertyView.Rows.Add(temp[i].Name, temp[i].GetValue(psm, null));
		//    }
		//    dataGridProperties.Items.Refresh();
		//    DrawPsm(psm.Ms2ScanNumber, psm.FullSequence);

		//    double wid = 0;
		//    dataGridProperties.HorizontalScrollBarVisibility = ScrollBarVisibility.Disabled;
		//    dataGridProperties.VerticalScrollBarVisibility = ScrollBarVisibility.Disabled;
		//    foreach (DataGridColumn col in dataGridProperties.Columns)
		//    {
		//        wid += col.ActualWidth;
		//    }
		//    PDFOutPut.Background = Brushes.White;
		//    PDFOutPut.ColumnDefinitions[0].Width = new GridLength(wid + 10);
		//    PDFOutPut.Measure(new Size(wid + gbPSM.ActualWidth + 10, 600));
		//    PDFOutPut.Arrange(new Rect(new Size(wid + gbPSM.ActualWidth + 10, 600)));
		//    dataGridProperties.Measure(new Size(wid + 22, 600));
		//    dataGridProperties.Arrange(new Rect(new Size(wid + 5, 600)));

		//    dataGridProperties.Arrange(new Rect(new Size(wid + 5, 600)));
		//    var rtb = new RenderTargetBitmap((int)(wid + gbPSM.ActualWidth) + 11, 600, 96, 96, PixelFormats.Pbgra32);

		//    rtb.Render(PDFOutPut);
		//    BitmapFrame bf = BitmapFrame.Create(rtb);

		//    var encoder = new BmpBitmapEncoder();
		//    encoder.Frames.Add(bf);
		//    using (var stream = new MemoryStream())
		//    {
		//        encoder.Save(stream);
		//        var img = System.Drawing.Image.FromStream(stream);
		//        PdfWriter.WriteToPdf(img, (int)(wid + gbPSM.ActualWidth) + 11, 600, path);
		//    }

		//    dataGridProperties.HorizontalScrollBarVisibility = ScrollBarVisibility.Visible;
		//    dataGridProperties.VerticalScrollBarVisibility = ScrollBarVisibility.Visible;
		//    PDFOutPut.ColumnDefinitions[0].Width = new GridLength(1, GridUnitType.Star);
		//}
	};
}
