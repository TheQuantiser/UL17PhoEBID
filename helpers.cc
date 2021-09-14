#include "map"
#include "string"
#include "TF1.h"
#include "iostream"
#include "fstream"
#include <dirent.h>
#include <sys/stat.h>


std::string removeNonAlpha(std::string word);
Char_t isDirectory(std::string filePath, Bool_t _verbose);
std::string getFileName(std::string _filepath);
Bool_t file_exists(std::string fileName);
std::vector<std::string> split_string(std::string _string, std::string _delimiter, Bool_t _trim);
void trim(std::string &s);
void ltrim(std::string &s);
void rtrim(std::string &s);


class CSVReader;
struct isoCorrMap;


std::string removeNonAlpha(std::string word) {
	word.erase(std::remove_if(word.begin(), word.end(),
	[](char ch) {
		return !::iswalnum(ch);
	}), word.end());
	return word;
};

Char_t isDirectory(std::string filePath, Bool_t _verbose=0) {
	DIR* dir = opendir(filePath.c_str());
	if (dir) {
		closedir(dir);
		return 1;
	} else if (ENOENT == errno) {
		return 0;
	} else {
		if (_verbose) std::cout<<"Error checking path : "<<filePath<<std::endl;
		return -1;
	}
}


std::string getFileName(std::string _filepath) {
	constexpr char sep = '/';
	size_t i = _filepath.rfind(sep, _filepath.length());
	if (i != string::npos) {
		return (_filepath.substr(i+1, _filepath.length() - i));
	}
	return (_filepath);
};



Bool_t file_exists(std::string fileName) {
	// std::ifstream infile(fileName);
	// return infile.good();
	if (isDirectory(fileName)==1) return 0;
	struct stat buffer;
	return (stat (fileName.c_str(), &buffer) == 0);
};


std::vector<std::string> split_string(std::string _string, std::string _delimiter=",", Bool_t _trim=1) {
	size_t pos = 0;
	std::string token;
	std::vector<std::string> res;
	while ((pos = _string.find(_delimiter)) != std::string::npos) {
		token = _string.substr(0, pos);
		if (_trim) trim(token);
		_string.erase(0, pos + _delimiter.length());
		res.push_back(token);
	}
	if (_trim) trim(_string);
	res.push_back(_string);
	return res;
};


void trim(std::string &s) {
	ltrim(s);
	rtrim(s);
}

void ltrim(std::string &s) {
	s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
		return !std::isspace(ch);
	}));
}



void rtrim(std::string &s) {
	s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
		return !std::isspace(ch);
	}).base(), s.end());
}


class CSVReader {
  private:

	void parseFile(std::string fileName, std::string delimeter = ",") {
		if (!file_exists(fileName)) {
			std::cout<<"Error! Cannot CSV-parse! File does not exist:"<<fileName<<std::endl;
		}
		dataList.clear();
		std::ifstream file(fileName);
		std::string line = "";
		while (getline(file, line)) {
			std::vector<std::string> line_split = split_string(line, delimeter);
			if (line_split.empty()) continue;
			dataList.push_back(line_split);
		}
		file.close();
	};

  public:
	CSVReader() {};

	CSVReader(std::string fileName, std::string delimeter = ",") {
		parseFile(fileName, delimeter);
	};

	std::vector<std::vector<std::string> > getData() {
		return dataList;
	};

	std::vector<std::vector<std::string> > dataList;
};


struct isoCorrMap {

	std::map<Float_t, TF1*> theMap;
	std::map<Float_t, Float_t> zeroVals;
	std::string mapFile;

	isoCorrMap() {};

	isoCorrMap(std::string mapFile, Int_t fColumn = 3, Bool_t _verbose=1, std::string _delimiter=";", Int_t _firstLine = 0) {
		init(mapFile, fColumn, _verbose, _delimiter, _firstLine);
	};

	~isoCorrMap() {
		for (std::map<Float_t, TF1*>::iterator iEl 		= 	theMap.begin(); iEl != theMap.end(); iEl++) {
			delete iEl->second;
			// iEl->second = nullptr;
		}

		theMap.clear();
	};

	void init(std::string _mapFile, Int_t fColumn = 3, Bool_t _verbose=1, std::string _delimiter=";", Int_t _firstLine = 0) {
		if (!file_exists(_mapFile)) {
			std::cout<<"\tError! Isolation map file "<< _mapFile<<" does not exist! "<<mapFile<<std::endl;
			exit(EXIT_FAILURE);
		}

		mapFile = _mapFile;
		CSVReader readFile(mapFile, _delimiter);
		std::vector<std::vector<std::string>> isCorrData 		= 	readFile.getData();

		std::cout<<"Loading isolation corrections from "<<_mapFile<<std::endl;

		for (UInt_t iRow = _firstLine; iRow < isCorrData.size(); iRow++) {
			Float_t uBound 										= 	std::stof(isCorrData[iRow][1]);

			if (_verbose) {
				std::cout<<"\t\t\t"<<uBound<<"\t\t"<<isCorrData[iRow][fColumn]<<std::endl;
			}

			TF1* isCorrFunc 									= 	new TF1(removeNonAlpha(getFileName(mapFile) + isCorrData[iRow][1]).c_str(), isCorrData[iRow][fColumn].c_str(), 0.,
			        13000.);
			theMap[uBound] 										= 	isCorrFunc;

			zeroVals[uBound]									=	isCorrFunc->Eval(0);
		}

	};

	Float_t getIsoCorr(Float_t absEta, Float_t sendaryEn, Bool_t _keepZero = 0) {
		if (theMap.empty()) return 0.;
		std::map<Float_t, TF1*>::iterator iBin = theMap.lower_bound(absEta);
		std::map<Float_t, Float_t>::iterator iBinZero = zeroVals.lower_bound(absEta);

		if (iBin == theMap.end()) {
			// std::cout<<"Error! Eta "<<absEta<<" is outside range specified in map "<<mapFile<<std::endl;
			return 0.;
		}
		if (_keepZero) return iBin->second->Eval(sendaryEn);
		return (iBin->second->Eval(sendaryEn) - iBinZero->second);
	};

	Float_t getEffectiveAreaAbs(Float_t eta, Float_t sendaryEn, Bool_t _keepZero = 0) {
		if (theMap.empty()) return 0.;
		Float_t 		absEta = std::abs(eta);
		std::map<Float_t, TF1*>::iterator iBin = theMap.lower_bound(absEta);
		std::map<Float_t, Float_t>::iterator iBinZero = zeroVals.lower_bound(absEta);

		if (iBin == theMap.end()) {
			std::cout<<"Error! Eta "<<eta<<" is outside range specified in map "<<mapFile<<std::endl;
			return 0.;
		}
		if (_keepZero) return iBin->second->Eval(sendaryEn);
		return (iBin->second->Eval(sendaryEn) - iBinZero->second);

	};
};
 
