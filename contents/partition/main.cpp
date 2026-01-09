#include <cstdio>
#include <iostream>
#include <string>
#include "oeis.h"
#include "debug.h"
#include "types.h"
#include "partition/alternate_equal.h"
#include "partition/condition_on_adjacent_diff.h"

int ri() {
	int n;
	scanf("%d", &n);
	return n;
}

int main(int argc, char **argv) {
	OEISGenerator gen;
	int seq_num = gen.get_sequence_number(argc, argv);
	if (seq_num == -1) {
		std::cerr << "No valid sequence number found" << std::endl;
		return 1;
	}
	OEISGenerator::GenOption opt = gen.get_cmdline_options(argc, argv, seq_num);
	
	opt.n_min = 1;
	opt.n_max = ri();
	opt.thread_num = 1;
	
	switch (seq_num) {
		// alternate equal
		case 351003 : gen.solve<A351003>(opt); break;
		case 351004 : gen.solve<A351004>(opt); break;
		case 351005 : gen.solve<A351005>(opt); break;
		case 351006 : gen.solve<A351006>(opt); break;
		case 351007 : gen.solve<A351007>(opt); break;
		case 351008 : gen.solve<A351008>(opt); break;
		case 351012 : gen.solve<A351012>(opt); break;
		// condition on adjacent diff
		case   4250 : gen.solve<A004250>(opt); break;
		case 120641 : gen.solve<A120641>(opt); break;
		case 323088 : gen.solve<A323088>(opt); break;
		case 323089 : gen.solve<A323089>(opt); break;
		case 323092 : gen.solve<A323092>(opt); break;
		case 323093 : gen.solve<A323093>(opt); break;
		case 323094 : gen.solve<A323094>(opt); break;
		case 350837 : gen.solve<A350837>(opt); break;
		case 350839 : gen.solve<A350839>(opt); break;
		case 350840 : gen.solve<A350840>(opt); break;
		case 350842 : gen.solve<A350842>(opt); break;
		case 350844 : gen.solve<A350844>(opt); break;
		case 350846 : gen.solve<A350846>(opt); break;
		default : std::cerr << "Unknown sequence : A" << seq_num << " found" << std::endl;
	}
	
	return 0;
}
