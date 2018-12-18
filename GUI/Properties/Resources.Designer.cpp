#include "Resources.Designer.h"

namespace MetaMorpheusGUI
{
	namespace Properties
	{
System::Resources::ResourceManager *Resources::resourceMan;
System::Globalization::CultureInfo *Resources::resourceCulture;

		Resources::Resources()
		{
		}

		System::Resources::ResourceManager *Resources::getResourceManager()
		{
			if (std::any::ReferenceEquals(resourceMan, nullptr))
			{
				System::Resources::ResourceManager *temp = new System::Resources::ResourceManager(L"MetaMorpheusGUI.Properties.Resources", Resources::typeid->Assembly);
				resourceMan = temp;

//C# TO C++ CONVERTER TODO TASK: A 'delete temp' statement was not added since temp was assigned to a field. Handle memory management manually.
			}
			return resourceMan;
		}

		System::Globalization::CultureInfo *Resources::getCulture()
		{
			return resourceCulture;
		}

		void Resources::setCulture(System::Globalization::CultureInfo *value)
		{
			resourceCulture = value;
		}
	}
}
