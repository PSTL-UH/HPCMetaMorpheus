#pragma once

#include <string>

using namespace iTextSharp::text;
using namespace iTextSharp::text::pdf;

namespace MetaMorpheusGUI
{
	class PdfWriter
	{
	public:
		static void WriteToPdf(System::Drawing::Image *img, double wid, double hei, const std::wstring &filePath);
	};
}
