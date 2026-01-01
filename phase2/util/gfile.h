#ifndef GINPUTFILE_H
#define GINPUTFILE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <string>
#include <fstream>

using namespace std;

#define MAX_LINE_LENGTH 32768

// Variable globale pour le débogage
class GLogFile ;

extern GLogFile glog ;


//------------------------------------------------------------
/*!
 *  \class GInputFile
 *
 *  \brief
 */

class GInputFile {
  private:
	std::string fileName;
    FILE* file;

    char line[MAX_LINE_LENGTH];
    char token[MAX_LINE_LENGTH];
    int indiceToken;
    int strlenline;

  public:
	GInputFile():file(NULL){}

    GInputFile(const char* nom):file(NULL){ fileName = std::string(nom) ;  }
	GInputFile(const std::string& nom):file(NULL){ fileName = nom ;  }

	inline bool isOpen(){return file!=NULL;}

    ~GInputFile(){}

    const char* getFileName(){return fileName.c_str();}
    void setFileName(char* st){fileName=st;}

    char* readLine() ;
    char* readUncommentedLine() ;
    char* readInt(int&) ;
    char* readIntInt(int &, int &) ;
    char* readIntIntInt(int &,int &,int &) ;
    char* readDoubleDoubleDouble(double &, double &, double &) ;

    char* getNextToken() ;
    char* getNextToken(char, bool canReturnEmpty=false) ;
    char getNextChar() ;

    string getEndOfLine();
    int   getNextIntToken() ;
    float getNextFloatToken() ;
    char* seek(char *) ;

    void open() ;
    void open(char *) ;
    void close() ;
};

//------------------------------------------------------------
/*!
 *  \class GOutputFile
 *
 * file that is used to write data. It has to be opened
 * before the data are written and closed after that.
 *
 */
class GOutputFile {
  protected :
    std::string fileName ;
    std::ofstream file ;
	int opened ;

  public:
	GOutputFile():opened(0){}
    GOutputFile(const std::string& nom):fileName(nom), opened(0){file<<std::fixed;}
	GOutputFile(char* nom):opened(0){fileName = std::string (nom) ; file<<std::fixed;}
	//GOutputFile(const char* nom):opened(0){fileName = std::string (nom) ;}
    virtual ~GOutputFile(){}

	virtual inline void setFileName (char* nom) { fileName = std::string (nom) ; }
	virtual inline void setFileName (const string& nom) { fileName = std::string (nom) ; }

	void open() {
		file.open (fileName.c_str(), std::ios_base::app|std::ios_base::out) ;
		opened = true ;
	}

	void openFromStart() {
		file.open (fileName.c_str(), std::ios_base::out) ;
		opened = true ;
	}

	void close() {
		if (opened) {
			file.close() ;
			opened = false ;
		}
	}

	virtual GOutputFile& operator<< (std::string& x) {file << x ; return *this ; }
	virtual GOutputFile& operator<< (char*  x) {file << x ; return *this ; }
	virtual GOutputFile& operator<< (int  x) {file << x ; return *this ; }
	virtual GOutputFile& operator<< (float x) {file << x ; return *this ; }
	virtual GOutputFile& operator<< (const char*  x) {file << x ; return *this ; }

} ;


//------------------------------------------------------------
/*!
 *  \class GLogFile
 *
 * file that is used to write data. It is automatically open/close
 * for each writing of data.
 * Advantage : all data are written even if there is a crash
 * Misadvantage : a lot of data writings takes much more time
 * than the GOutputFile
 *
 */
class GLogFile : public GOutputFile {
  public:
	GLogFile(){}
    GLogFile(const std::string& nom):GOutputFile(nom){}
	GLogFile(char* nom):GOutputFile(nom){}

	virtual GLogFile& operator<< (const std::string& x) {
		std::ofstream f (fileName.c_str(), std::ios_base::app|std::ios_base::out) ;
        f << x ;
        f.close();
		return *this ;
	}
	virtual GLogFile& operator<<(int x) {
		std::ofstream f (fileName.c_str(), std::ios_base::app|std::ios_base::out) ;
        f << x ;
        f.close();
		return *this ;
	}
	virtual GLogFile& operator<<(float x) {
		std::ofstream f (fileName.c_str(), std::ios_base::app|std::ios_base::out) ;
        f << x ;
        f.close();
		return *this ;
	}
	virtual GLogFile& operator<<(double x) {
		std::ofstream f (fileName.c_str(), std::ios_base::app|std::ios_base::out) ;
        f << x ;
        f.close();
		return *this ;
	}
	virtual GLogFile& operator<<(char* x) {
		std::ofstream f (fileName.c_str(), std::ios_base::app|std::ios_base::out) ;
        f << x ;
        f.close();
		return *this ;
	}

	virtual GLogFile& debug (const std::string& x) {
		std::ofstream f (fileName.c_str(), std::ios_base::app|std::ios_base::out) ;
		f << "Debug : " << x ;
        f.close();
		return *this ;
	}
	inline virtual GLogFile& debug (const char* x) {return debug(string(x)) ;}

	virtual GLogFile& info (const std::string& x) {
		std::ofstream f (fileName.c_str(), std::ios_base::app|std::ios_base::out) ;
		f << "Info : " << x ;
        f.close();
		return *this ;
	}
	inline virtual GLogFile& info (const char* x) {return info(string(x)) ;}

	virtual GLogFile& warning (const std::string& x) {
		std::ofstream f (fileName.c_str(), std::ios_base::app|std::ios_base::out) ;
		f << "Warning : " << x ;
        f.close();
		return *this ;
	}
	inline virtual GLogFile& warning (const char* x) {return warning(string(x)) ;}

	virtual GLogFile& error (const std::string& x) {
		std::ofstream f (fileName.c_str(), std::ios_base::app|std::ios_base::out) ;
		f << "ERROR : " << x ;
        f.close();
		return *this ;
	}
	inline virtual GLogFile& error (const char* x) {return error(string(x)) ;}

} ;


#endif
