#ifndef PARSE_COMMAND_H
#define PARSE_COMMAND_H

#include <string>
using std::string;

// ------------------------------------------------------------------------------------
// Parse int
int parseCommand( const string & line, const string & command, int & value, bool echo=true )
{
    int len = command.length();
    if( line.substr(0,len)==command )
    {
        sscanf(line.substr(len).c_str(),"%d",&value);
        if( echo ) printf("parseCommand: SETTING %s%d\n",command.c_str(),value);
        return 1;
    }
    return 0;
}

// Parse double
int parseCommand( const string & line, const string & command, double & value, bool echo=true )
{
    int len = command.length();
    if( line.substr(0,len)==command )
    {
        sscanf(line.substr(len).c_str(),"%le",&value);
        if( echo ) printf("parseCommand: SETTING %s%e\n",command.c_str(),value);
        return 1;
    }
    return 0;
}

// Parse string
int parseCommand( const string & line, const string & command, string & value, bool echo=true )
{
    int len = command.length();
    if( line.substr(0,len)==command )
    {
        value = line.substr(len);
        if( echo ) printf("parseCommand: SETTING %s%s\n",command.c_str(),value.c_str());
        return 1;
    }
    return 0;
}

#endif
