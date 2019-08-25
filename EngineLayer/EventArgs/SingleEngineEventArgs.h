#pragma once

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class MetaMorpheusEngine; }


namespace EngineLayer
{
	class SingleEngineEventArgs : public EventArgs
	{
	private:
		MetaMorpheusEngine *privateMyEngine;

	public:
		SingleEngineEventArgs(MetaMorpheusEngine *myEngine);

		MetaMorpheusEngine *getMyEngine() const;
		void setMyEngine(MetaMorpheusEngine *value);
	};
}
