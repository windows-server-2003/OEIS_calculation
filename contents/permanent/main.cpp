#define _GLIBCXX_DEBUG
#include <cstdio>
#include <iostream>
#include <string>
#include "oeis.h"
#include "debug.h"
#include "types.h"
#include "permanent/permanent_oeis.h"

int ri() {
	int n;
	scanf("%d", &n);
	return n;
}


int main(int argc, char **argv) {
	OEISGenerator gen;
	int seq_num = gen.get_sequence_number(argc, argv);
	OEISGenerator::GenOption opt = gen.get_cmdline_options(argc, argv, seq_num);
	
	if (seq_num == -1) {
		std::cerr << "No valid sequence number found" << std::endl;
		std::cerr << "Usage(for example): " << argv[0] << " A085719" << std::endl;
		return 1;
	}
	
	opt.n_min = 1;
	opt.n_max = ri();
	opt.thread_num = 4;
	opt.thread_limits = {{30, 2}, {35, 1}};
	
	switch (seq_num) {
		case  85719 : gen.solve<A085719>(opt); break;
		case  86759 : gen.solve<A086759>(opt); break;
		case  85807 : gen.solve<A085807>(opt); break;
		case 307783 : gen.solve<A307783>(opt); break;
		case 322277 : gen.solve<A322277>(opt); break;
		case 278847 : gen.solve<A278847>(opt); break;
		case 278925 : gen.solve<A278925>(opt); break;
		case 278926 : gen.solve<A278926>(opt); break;
		case 278927 : gen.solve<A278927>(opt); break;
		case 346934 : gen.solve<A346934>(opt); break;
		case 347768 : gen.solve<A347768>(opt); break;
		case 278845 : gen.solve<A278845>(opt); break;
		case 278857 : gen.solve<A278857>(opt); break;
		case 278858 : gen.solve<A278858>(opt); break;
		case 203264 : gen.solve<A203264>(opt); break;
		case 303000 : gen.solve<A303000>(opt); break;
		case 303001 : gen.solve<A303001>(opt); break;
		case 179079 : gen.solve<A179079>(opt); break;
		case 204234 : gen.solve<A204234>(opt); break;
		case 204235 : gen.solve<A204235>(opt); break;
		case 204236 : gen.solve<A204236>(opt); break;
		case 204239 : gen.solve<A204239>(opt); break;
		case 204241 : gen.solve<A204241>(opt); break;
		case 204248 : gen.solve<A204248>(opt); break;
		case 204249 : gen.solve<A204249>(opt); break;
		case 204251 : gen.solve<A204251>(opt); break;
		case 204252 : gen.solve<A204252>(opt); break;
		case 204254 : gen.solve<A204254>(opt); break;
		case 204256 : gen.solve<A204256>(opt); break;
		case 204258 : gen.solve<A204258>(opt); break;
		case 204262 : gen.solve<A204262>(opt); break;
		case 204264 : gen.solve<A204264>(opt); break;
		case 204265 : gen.solve<A204265>(opt); break;
		case 204268 : gen.solve<A204268>(opt); break;
		case 204422 : gen.solve<A204422>(opt); break;
		case 204424 : gen.solve<A204424>(opt); break;
		case 204426 : gen.solve<A204426>(opt); break;
		case 204428 : gen.solve<A204428>(opt); break;
		case 204430 : gen.solve<A204430>(opt); break;
		case 204432 : gen.solve<A204432>(opt); break;
		case 204434 : gen.solve<A204434>(opt); break;
		case 204436 : gen.solve<A204436>(opt); break;
		case 204438 : gen.solve<A204438>(opt); break;
		case 204440 : gen.solve<A204440>(opt); break;
		case 204442 : gen.solve<A204442>(opt); break;
		case 204444 : gen.solve<A204444>(opt); break;
		case 204446 : gen.solve<A204446>(opt); break;
		case 204448 : gen.solve<A204448>(opt); break;
		case 204546 : gen.solve<A204546>(opt); break;
		case 204548 : gen.solve<A204548>(opt); break;
		case 204550 : gen.solve<A204550>(opt); break;
		case 322909 : gen.solve<A322909>(opt); break;
		case 323255 : gen.solve<A323255>(opt); break;
		case 330087 : gen.solve<A330087>(opt); break;
	}
	
	return 0;
}
