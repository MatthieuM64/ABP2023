/*LIBRARY FOR SOME SPECIAL FUNCTIONS*/

#ifndef DEF_SPECIAL_FUNCTIONS_MANGEAT_CPP
#define DEF_SPECIAL_FUNCTIONS_MANGEAT_CPP

using namespace std;

class timeProgram
{
	time_t begin_time;
	
	public:
	
	timeProgram();
	unsigned long int Time() const;
	unsigned long int BeginTime() const;
	string TimeRun(const string &c) const;
	string TimeRun() const;
};

timeProgram::timeProgram()
{
	time(&begin_time);
}

unsigned long int timeProgram::BeginTime() const
{
	struct tm y2k = {0};
	y2k.tm_hour = 0;   y2k.tm_min = 0; y2k.tm_sec = 0;
	y2k.tm_year = 115; y2k.tm_mon = 0; y2k.tm_mday = 1;
	unsigned long int time_from_2015=floor(difftime(begin_time,mktime(&y2k)));
	return time_from_2015;
}

unsigned long int timeProgram::Time() const
{
	time_t actual_time;
	time(&actual_time);
	unsigned long int time_sec=floor(difftime(actual_time,begin_time));
	return time_sec;
}

string timeProgram::TimeRun(const string &c) const									//Return the execution time
{
	char time[100];
	unsigned long int time_sec=Time();
	if (time_sec<60)
	{
		sprintf(time,"%s-ctime=00h00m%02ds",c.c_str(),int(time_sec));
	}
	else if (time_sec<3600)
	{
		sprintf(time,"%s-ctime=00h%02dm%02ds",c.c_str(),int(time_sec/60),int(time_sec%60));
	}
	else if (time_sec<86400)
	{
		sprintf(time,"%s-ctime=%02dh%02dm%02ds",c.c_str(),int(time_sec/3600),int((time_sec/60)%60),int(time_sec%60));
	}
	else
	{
		sprintf(time,"%s-ctime=%dd%02dh%02dm%02ds",c.c_str(),int(time_sec/86400),int((time_sec/3600)%24),int((time_sec/60)%60),int(time_sec%60));
	}
	string time_str(time);
	return time_str;
}

string timeProgram::TimeRun() const											//Return the execution time
{
	return TimeRun("");
}

timeProgram running_time;

double square(const double &r)
{
	return r*r;
}

double cube(const double &r)
{
	return r*r*r;
}

#endif
