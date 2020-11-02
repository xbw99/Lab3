//////////////////////////////////////////////////////////////////////
// C++ source file generated by SpecC V2.2.1
// Design: main
// File:   main.cc
// Time:   Mon Nov  2 16:44:43 2020
//////////////////////////////////////////////////////////////////////

// Note: User-defined include files are inlined in this file.

// Note: System-defined include files are inlined in this file.

// Note: Run-time debugging is enabled in this file.

#include "main.h"


unsigned int _IDcnt = 0;
// channel class definitions /////////////////////////////////////////

c_image_queue::c_image_queue(const unsigned long int (&size),  const char *_InstanceName, _specc::class_type *_Parent)
    : _specc::channel("c_image_queue", _InstanceName, _Parent), size(size),
    buffer(0),
    n(0ul),
    p(0ul),
    wr(0ul),
    ws(0ul)
{   
}

c_image_queue::~c_image_queue(void)
{   
}

#line 10 "./image_queue.sc"
void c_image_queue::cleanup(void) {const char *_m; _SetCurrentMethod("cleanup", &_m); if ( !n) { free(buffer); buffer = 0;
    }_RestoreCurrentMethod(_m);
}

#line 10 "./image_queue.sc"
void c_image_queue::receive(unsigned char (*d)[76800]) {_specc::class_type *_p; const char *_m; _SetCurrentInst("receive", &_p, &_m); while( !n) { wr++ ; _specc::wait(event(&r), ((void*)0)); wr-- ;
    }

#line 10 "./image_queue.sc"
    if (n <= p) { { unsigned int _scc_index_0; for(_scc_index_0=0;_scc_index_0<76800;_scc_index_0++) ( *d)[_scc_index_0] = (buffer[p - n])[_scc_index_0]; }
    }
    else 

#line 10 "./image_queue.sc"
    {    { unsigned int _scc_index_0; for(_scc_index_0=0;_scc_index_0<76800;_scc_index_0++) ( *d)[_scc_index_0] = (buffer[p + size - n])[_scc_index_0]; }
    }

#line 10 "./image_queue.sc"
    n-- ; if (ws) { _specc::notify(event(&s), ((void*)0));
    }

#line 10 "./image_queue.sc"
    cleanup();_RestoreCurrentInst(_p, _m);
}

#line 10 "./image_queue.sc"
void c_image_queue::send(unsigned char d[76800]) {_specc::class_type *_p; const char *_m; _SetCurrentInst("send", &_p, &_m); while(n >= size) { ws++ ; _specc::wait(event(&s), ((void*)0)); ws-- ;
    }

#line 10 "./image_queue.sc"
    setup(); { unsigned int _scc_index_0; for(_scc_index_0=0;_scc_index_0<76800;_scc_index_0++) (buffer[p])[_scc_index_0] = (d)[_scc_index_0]; } p++ ; if (p >= size) { p = 0;
    }

#line 10 "./image_queue.sc"
    n++ ; if (wr) { _specc::notify(event(&r), ((void*)0));
    }_RestoreCurrentInst(_p, _m);
}

#line 10 "./image_queue.sc"
void c_image_queue::setup(void) {const char *_m; _SetCurrentMethod("setup", &_m); if ( !buffer) { unsigned char dummy[76800]; unsigned long int i; if ( !(buffer = (unsigned char (*)[76800])malloc(sizeof(unsigned char [76800]) * size))) { perror("c_typed_queue"); abort();
	}

#line 10 "./image_queue.sc"
	for(i = 0; i < size; i++ ) { memcpy( &buffer[i],  &dummy, sizeof(unsigned char [76800]));
	}
    }_RestoreCurrentMethod(_m);
}

// behavior class definitions ////////////////////////////////////////

#line 86 "main.cc"
Monitor::Monitor(unsigned int _idcnt, i_image_receiver (&ch),  const char *_InstanceName, _specc::class_type *_Parent)
    : _specc::behavior(_idcnt, "Monitor" , _InstanceName, _Parent), ch(ch)
{   
}

Monitor::~Monitor(void)
{   
}

#line 14 "./monitor.sc"
void Monitor::main(void) {_specc::class_type *_p; const char *_m; _SetActiveInst("main", &_p, &_m);
    ch.receive( &IMAGE);

    outfilename = (char *)"ba.pgm";
    comment = (char *)"";
    maxval = 255;
    image = IMAGE;

    if (outfilename == ((void *)0)) fp = stdout;
    else  {
	if ((fp = fopen(outfilename, "w")) == ((void *)0)) {
	    fprintf(stderr, "Error writing the file %s in write_pgm_image().\n", 
		outfilename);
	    { _RestoreActiveInst(_p, _m);return ; }
	}
    }



    fprintf(fp, "P5\n%d %d\n", cols, rows);
    if (comment != ((void *)0))
	if (strlen(comment) <= 70) fprintf(fp, "# %s\n", comment);
    fprintf(fp, "%d\n", maxval);


    if ((unsigned int)rows != fwrite(image, cols, rows, fp)) {
	fprintf(stderr, "Error writing the image data in write_pgm_image().\n");
	if (fp != stdout) fclose(fp);
	{ _RestoreActiveInst(_p, _m);return ; }
    }

    if (fp != stdout) fclose(fp);
    { _RestoreActiveInst(_p, _m);return ; }
}

#line 132 "main.cc"
DataIn::DataIn(unsigned int _idcnt, i_image_receiver (&ch1), i_image_sender (&ch2),  const char *_InstanceName, _specc::class_type *_Parent)
    : _specc::behavior(_idcnt, "DataIn" , _InstanceName, _Parent), ch1(ch1), ch2(ch2)
{   
}

DataIn::~DataIn(void)
{   
}

#line 17 "./platform.sc"
void DataIn::main(void) {_specc::class_type *_p; const char *_m; _SetActiveInst("main", &_p, &_m);
    while(true) {
	ch1.receive( &IMAGE);
	ch2.send(IMAGE);
    }_RestoreActiveInst(_p, _m);
}

#line 150 "main.cc"
DataOut::DataOut(unsigned int _idcnt, i_image_receiver (&ch1), i_image_sender (&ch2),  const char *_InstanceName, _specc::class_type *_Parent)
    : _specc::behavior(_idcnt, "DataOut" , _InstanceName, _Parent), ch1(ch1), ch2(ch2),
    counter(0)
{   
}

DataOut::~DataOut(void)
{   
}

#line 29 "./platform.sc"
void DataOut::main(void) {_specc::class_type *_p; const char *_m; _SetActiveInst("main", &_p, &_m);
    int stop = 10;
    while(true) {
	ch1.receive( &IMAGE);
	ch2.send(IMAGE);
	if (counter++  == stop) {
	    exit(1);
	}
    }_RestoreActiveInst(_p, _m);
}

#line 173 "main.cc"
DUT::DUT(unsigned int _idcnt, i_image_receiver (&ch1), i_image_sender (&ch2),  const char *_InstanceName, _specc::class_type *_Parent)
    : _specc::behavior(_idcnt, "DUT" , _InstanceName, _Parent), ch1(ch1), ch2(ch2)
{   
}

DUT::~DUT(void)
{   
}

#line 43 "./platform.sc"
void DUT::main(void) {_specc::class_type *_p; const char *_m; _SetActiveInst("main", &_p, &_m);
    ch1.receive( &IMAGE);
    ch2.send(IMAGE);_RestoreActiveInst(_p, _m);
}

#line 189 "main.cc"
Platform::Platform(unsigned int _idcnt, i_image_receiver (&toPlatform), i_image_sender (&outPlatform),  const char *_InstanceName, _specc::class_type *_Parent)
    : _specc::behavior(_idcnt, "Platform" , _InstanceName, _Parent), toPlatform(toPlatform), outPlatform(outPlatform),
    _scc_const_port_0(2ul),
    _scc_const_port_1(2ul),
    canny(_IDcnt, inDut, outDut, "canny", this),
    din(_IDcnt, toPlatform, inDut, "din", this),
    dout(_IDcnt, outDut, outPlatform, "dout", this),
    inDut(_scc_const_port_0, "inDut", this),
    outDut(_scc_const_port_1, "outDut", this)
{   
}

Platform::~Platform(void)
{   
}

#line 56 "./platform.sc"
void Platform::main(void) {_specc::class_type *_p; const char *_m; _SetActiveInst("main", &_p, &_m);
    { _specc::fork _scc_fork_0(&din), _scc_fork_1(&canny), _scc_fork_2(&dout); _specc::par(&_scc_fork_0, &_scc_fork_1, &_scc_fork_2, ((_specc::fork*)0));
    }

#line 57 "./platform.sc"
    ;_RestoreActiveInst(_p, _m);
}

#line 215 "main.cc"
Stimulus::Stimulus(unsigned int _idcnt, i_image_sender (&ch), char *(&infilename),  const char *_InstanceName, _specc::class_type *_Parent)
    : _specc::behavior(_idcnt, "Stimulus" , _InstanceName, _Parent), ch(ch), infilename(infilename),
    image(0)
{   
}

Stimulus::~Stimulus(void)
{   
}

#line 11 "./stimulus.sc"
void Stimulus::main(void) {_specc::class_type *_p; const char *_m; _SetActiveInst("main", &_p, &_m);

    struct _IO_FILE *fp;
    char buf[71];


    if (infilename == ((void *)0)) fp = stdin;
    else  {
	if ((fp = fopen(infilename, "r")) == ((void *)0)) {
	    fprintf(stderr, "Error reading the file %s in read_pgm_image().\n", 
		infilename);
	    { _RestoreActiveInst(_p, _m);return ; }
	}
    }

    fgets(buf, 70, fp);
    if (strncmp(buf, "P5", 2) != 0) {
	fprintf(stderr, "The file %s is not in PGM format in ", infilename);
	fprintf(stderr, "read_pgm_image().\n");
	if (fp != stdin) fclose(fp);
	{ _RestoreActiveInst(_p, _m);return ; }
    }

    do  { fgets(buf, 70, fp);
    }
    while(

#line 34 "./stimulus.sc"
	buf[0] == '#');
    __isoc99_sscanf(buf, "%d %d",  &cols,  &rows);
    do  { fgets(buf, 70, fp);
    }
    while(

#line 36 "./stimulus.sc"
	buf[0] == '#');


    image = IMAGE;

    if ((unsigned int)rows != fread(( *( &image)), cols, rows, fp)) {
	fprintf(stderr, "Error reading the image data in read_pgm_image().\n");
	if (fp != stdin) fclose(fp);
	{ _RestoreActiveInst(_p, _m);return ; }
    }

    if (fp != stdin) fclose(fp);

    ch.send(IMAGE);
    { _RestoreActiveInst(_p, _m);return ; }
}

#line 279 "main.cc"
Main::Main(unsigned int _idcnt,  const char *_InstanceName, _specc::class_type *_Parent)
    : _specc::class_type(_idcnt, "Main" , _InstanceName, _Parent),
    _scc_const_port_0(2ul),
    _scc_const_port_1(2ul),
    mon(_IDcnt, outPlatform, "mon", this),
    outPlatform(_scc_const_port_0, "outPlatform", this),
    plat(_IDcnt, toPlatform, outPlatform, "plat", this),
    stim(_IDcnt, toPlatform, infilename, "stim", this),
    toPlatform(_scc_const_port_1, "toPlatform", this)
{   
}

Main::~Main(void)
{   
}

#line 14 "main.sc"
int Main::main(int argc, char *argv[]) {_specc::class_type *_p; const char *_m; _SetActiveInst("main", &_p, &_m);
    infilename = argv[1];

    { _specc::fork _scc_fork_0(&stim), _scc_fork_1(&plat), _scc_fork_2(&mon); _specc::par(&_scc_fork_0, &_scc_fork_1, &_scc_fork_2, ((_specc::fork*)0));
    }

#line 17 "main.sc"
    ;
    { _RestoreActiveInst(_p, _m);return 0; }
}

#line 308 "main.cc"
Main _scc_main(_IDcnt,  "Main", 0);

int main(int argc, char *argv[])
{   
    int _scc_main_return;
    
    _specc::start();
    _scc_main_return = _scc_main.main(argc, argv);
    _specc::end();
    return(_scc_main_return);
}

void _scc_bit4_err_handle(
    const _bit4& bit4vec)
{   
    char temp_bits[1024], *p;
    p=bit2str(2,&temp_bits[1023], bit4vec);
    _specc::abort(
	"ERROR:\t Casting a bit4 vector failed \n"
	"Bit4 vector contains X/Z values %s\n"
	"Simulation aborted.\n", p);
	
}

//////////////////////////////////////////////////////////////////////
// End of file main.cc
//////////////////////////////////////////////////////////////////////