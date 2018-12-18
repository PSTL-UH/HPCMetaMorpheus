#pragma once

#include <string>


namespace MetaMorpheusGUI
{
	class BaseDraw
	{
		/// <summary>
		/// Draw the line seperator @ top
		/// </summary>
		/// <param name="cav">Canvas to where draws the line</param>
		/// <param name="topLoc">Location of the starting point on top (exp: 0,0)</param>
	public:
		static void topSplittingDrawing(Canvas *cav, Point *topLoc, Color *clr, const std::wstring &footnote);

		/// <summary>
		/// Draw the line seperator @ bottom
		/// </summary>
		/// <param name="cav">Canvas to where draws the line</param>
		/// <param name="botLoc">Location of the starting point on bottom (exp: 0,50)</param>
		static void botSplittingDrawing(Canvas *cav, Point *botLoc, Color *clr, const std::wstring &footnote);

		/// <summary>
		/// Create text blocks on canvas
		/// </summary>
		/// <param name="cav"> Canvas Board </param>
		/// <param name="loc"> Provate the (x,y) coordinates for textblock</param>
		/// <param name="txt"> Message for textblock</param>
		/// <returns> the width of current addup</returns>
		static void txtDrawing(Canvas *cav, Point *loc, const std::wstring &txt, Brush *clr);

		/// <summary>
		/// Display texts surounded by circle
		/// </summary>
		/// <param name="cav"></param>
		/// <param name="loc"></param>
		/// <param name="txt"></param>
		/// <param name="clr"></param>
		/// <returns></returns>
		static void circledTxtDraw(Canvas *cav, Point *loc, SolidColorBrush *clr);

		/// <summary>
		/// Clear canvas board
		/// </summary>
		/// <param name="cav">board to clear</param>
		static void clearCanvas(Canvas *cav);
	};
}
