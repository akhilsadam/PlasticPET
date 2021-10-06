#define PBWIDTH 20
#include <stdio.h>
#include <stdlib.h>
#include "G4Run.hh"

class RunTools
{
private:
	static int currentEventNumber;
	static int totalEvents;
public:
	RunTools(const G4Run* run) {
		currentEventNumber = 0;
		totalEvents = (int)run->GetNumberOfEventToBeProcessed();
	}
	static void completedEvent() {
		currentEventNumber = currentEventNumber + 1;
		printProgress();
	}
	static int GetTotalEvents() {
		return totalEvents;
	}
private:
	static void printProgress() {
		double percentage = (double)(currentEventNumber / totalEvents);
		int lpad = (int)(percentage * PBWIDTH);
		int rpad = PBWIDTH - lpad;
		//char b = 219;
		// system("COLOR 09");
		//printf("\r%3d /%3d [%.*s%*s]", currentEventNumber, totalEvents, lpad, PBSTR, rpad, "");
		printf("\r Events:%3d /%3d [", currentEventNumber, totalEvents);
		for (int i = 0; i < lpad; i++) {
			//printf("%c", b);
			printf("-");
		}
		for (int i = 0; i < rpad; i++) {
			printf(" ");
		}
		printf("]");
		fflush(stdout);
	}
};
