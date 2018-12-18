#include "BaseSequenceAnnotation.h"


namespace MetaMorpheusGUI
{

	void BaseDraw::topSplittingDrawing(Canvas *cav, Point *topLoc, Color *clr, const std::wstring &footnote)
	{
		double x = topLoc->X, y = topLoc->Y;
		Polyline *bot = new Polyline();
		bot->Points = new PointCollection()
		{
			new Point(x + 10, y),
			new Point(x, y + 10),
			new Point(x, y + 40)
		};
		bot->Stroke = new SolidColorBrush(clr);
		bot->StrokeThickness = 2;
		cav->Children->Add(bot);
		Canvas::SetZIndex(bot, 1); //on top of any other things in canvas
		TextBlock *tb = new TextBlock();
		tb->Foreground = new SolidColorBrush(clr);
		tb->Text = footnote;
		tb->FontSize = 10;
		Canvas::SetTop(tb, y - 10);
		Canvas::SetLeft(tb, x + 10);
		Canvas::SetZIndex(tb, 2);
		cav->Children->Add(tb);

//C# TO C++ CONVERTER TODO TASK: A 'delete tb' statement was not added since tb was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete bot' statement was not added since bot was passed to a method or constructor. Handle memory management manually.
	}

	void BaseDraw::botSplittingDrawing(Canvas *cav, Point *botLoc, Color *clr, const std::wstring &footnote)
	{
		double x = botLoc->X, y = botLoc->Y;
		Polyline *bot = new Polyline();
		bot->Points = new PointCollection()
		{
			new Point(x - 10, y),
			new Point(x, y - 10),
			new Point(x, y - 40)
		};
		bot->Stroke = new SolidColorBrush(clr);
		bot->StrokeThickness = 2;
		Canvas::SetZIndex(bot, 1); //on top of any other things in canvas
		cav->Children->Add(bot);
		TextBlock *tb = new TextBlock();
		tb->Foreground = new SolidColorBrush(clr);
		tb->Text = footnote;
		tb->FontSize = 10;
		Canvas::SetTop(tb, y - 8);
		Canvas::SetLeft(tb, x - 22);
		Canvas::SetZIndex(tb, 2);
		cav->Children->Add(tb);

//C# TO C++ CONVERTER TODO TASK: A 'delete tb' statement was not added since tb was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete bot' statement was not added since bot was passed to a method or constructor. Handle memory management manually.
	}

	void BaseDraw::txtDrawing(Canvas *cav, Point *loc, const std::wstring &txt, Brush *clr)
	{
		TextBlock *tb = new TextBlock();
		tb->Foreground = clr;
		tb->Text = txt;
		tb->Height = 30;
		tb->FontSize = 25;
		tb->FontWeight = FontWeights::Bold;
		tb->FontFamily = new FontFamily(L"Courier New"); // monospaced font

		Canvas::SetTop(tb, loc->Y);
		Canvas::SetLeft(tb, loc->X);
		Panel::SetZIndex(tb, 2); //lower priority
		cav->Children->Add(tb);
		cav->UpdateLayout();

//C# TO C++ CONVERTER TODO TASK: A 'delete tb' statement was not added since tb was passed to a method or constructor. Handle memory management manually.
	}

	void BaseDraw::circledTxtDraw(Canvas *cav, Point *loc, SolidColorBrush *clr)
	{
		Ellipse *circle = new Ellipse();
		circle->Width = 24;
		circle->Height = 24;
		circle->Stroke = clr;
		circle->StrokeThickness = 1;
		circle->Fill = clr;
		circle->Opacity = 0.7;
		Canvas::SetLeft(circle, loc->X);
		Canvas::SetTop(circle, loc->Y);
		Panel::SetZIndex(circle, 1);
		cav->Children->Add(circle);

//C# TO C++ CONVERTER TODO TASK: A 'delete circle' statement was not added since circle was passed to a method or constructor. Handle memory management manually.
	}

	void BaseDraw::clearCanvas(Canvas *cav)
	{
		cav->Children->Clear();
	}
}
