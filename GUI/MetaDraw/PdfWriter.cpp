#include "PdfWriter.h"

using namespace iTextSharp::text;
using namespace iTextSharp::text::pdf;

namespace MetaMorpheusGUI
{

	void PdfWriter::WriteToPdf(System::Drawing::Image *img, double wid, double hei, const std::wstring &filePath)
	{
		Rectangle tempVar(static_cast<float>(wid), static_cast<float>(hei));
		Document *doc = new Document(&tempVar);
		FileStream tempVar2(filePath, FileMode::Create);
		iTextSharp::text::pdf::PdfWriter::GetInstance(doc, &tempVar2);
		doc->Open();
		Image *pdfImage = Image::GetInstance(img, System::Drawing::Imaging::ImageFormat::Bmp);
		pdfImage->SetAbsolutePosition(0, 0);
		doc->Add(pdfImage);
		doc->Close();

//C# TO C++ CONVERTER TODO TASK: A 'delete doc' statement was not added since doc was passed to a method or constructor. Handle memory management manually.
	}
}
