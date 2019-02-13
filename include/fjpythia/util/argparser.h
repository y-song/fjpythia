#ifndef FJPY_ARGPARSER_HH
#define FJPY_ARGPARSER_HH

// modified from https://stackoverflow.com/questions/865668/how-to-parse-command-line-arguments-in-c

#include <string>
#include <vector>

namespace FJPyUtil
{
	class ArgParser
	{
		public:
			static ArgParser& Instance()
			{
				static ArgParser _instance; // Guaranteed to be destroyed. // Instantiated on first use.
				return _instance;
			}

			void addOpts(int &argc, char **argv)
			{
				for (int i=1; i < argc; ++i)
					this->tokens.push_back(std::string(argv[i]));
			}

			virtual ~ArgParser();

			const std::string& getOpt(const std::string &option)
			{
				registerOpt(option);
				std::vector<std::string>::const_iterator itr;
				itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
				if (itr != this->tokens.end() && ++itr != this->tokens.end())
				{
					return *itr;
				}
				static const std::string empty_string("");
				return empty_string;
			}

			const std::string& getOpt(const char* option)
			{
				std::string opt(option);
				return getOpt(opt);
			}

			bool isSet(const std::string &option)
			{
				registerOpt(option);
				return std::find(this->tokens.begin(), this->tokens.end(), option)
					   != this->tokens.end();
			}

			bool isSet(const char* option)
			{
				std::string opt(option);
				return isSet(opt);
			}

		private:
			ArgParser(); // private construct
			ArgParser(const ArgParser&) = delete; // non-copyable
			ArgParser(ArgParser&&) = delete; // and non-movable
			// as there is only one object, assignment would always be assign to self
			ArgParser& operator= ( const ArgParser& ) = delete ;
			ArgParser& operator= ( ArgParser&& ) = delete ;

			void registerOpt(const std::string &option)
			{
				if (std::find(this->regopts.begin(), this->regopts.end(), option) == this->regopts.end())
					regopts.push_back(option);
			}

			std::vector <std::string> tokens;
			std::vector <std::string> regopts; // registered/defined/querried in the code opts
	};
}

#endif