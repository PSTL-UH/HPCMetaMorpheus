#pragma once

#include <string>
#include <stdexcept>
#include <any>
#include "stringbuilder.h"
#include "tangible_filesystem.h"

using namespace EngineLayer;
using namespace Newtonsoft::Json::Linq;

namespace MetaMorpheusGUI
{
	/// <summary>
	/// Interaction logic for MetaUpdater.xaml
	/// </summary>
	class MetaUpdater : public Window
	{
	public:
		MetaUpdater();

//C# TO C++ CONVERTER TODO TASK: Methods returning tuples are not converted by C# to C++ Converter:
//		public (int, int, int) GetVersionNumber(string VersionNode)
//		{
//			try
//			{
//				var split = VersionNode.Split('.');
//
//				return (int.Parse(split[0]), int.Parse(split[1]), int.Parse(split[2]));
//			}
//			catch (FormatException)
//			{
//				return (0, 0, 0);
//			}
//		}

	private:
		void InstallerClicked(std::any sender, RoutedEventArgs *e);

		void ReleaseHandler();

		void PortableClicked(std::any semder, RoutedEventArgs *e);

		void NoClicked(std::any semder, RoutedEventArgs *e);
	};
}
