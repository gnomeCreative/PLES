#ifndef MACROS_H_
#define MACROS_H_

// ATTENTION!!! MAJOR MODIFICATION HERE!
// instead of cout here we had *os


// compile options
//#define MEASURE_PERFORMANCE
#define SHELL_OMP

// a useful print info macro
#define INFO __FILE__ << ":" << __LINE__ << " "

// assert a condition
#define ASSERT(cond) if(!(cond)) { cout << INFO << "ERROR: Condition " #cond " not fulfilled!" << std::endl; exit(1); }
#define REQUEST(cond) if(!(cond)) { cout << INFO << "WARNING: Condition " #cond " not fulfilled!" << std::endl; }

// read a parameter from the config file, and if it is also given on the command line, overwrite it with that one
#define BASE_PARSE_CFG(cfg,type,var,name,def,def_typed) \
		if (!cfg.have_variable(name) && !command_line.search("-" name))\
		cout << INFO << "WARNING: Parameter " name " is missing! Using default: " << def << std::endl;\
		type var = cfg(name, def_typed);\
		if (command_line.search("-" name))\
		var = command_line.next(var);

//// read a 3-vector parameter from the config file, and if it is also given on the command line, overwrite it with that one
//#define BASE_PARSE_CFG(cfg,type,var,name,def,def_typed) \
//  if (!cfg.have_variable(name) && !command_line.search("-" name))\
//    cout << INFO << "WARNING: Parameter " name " is missing! Using default: " << def << std::endl;\
//  type var = cfg(name, def_typed);\
//  if (command_line.search("-" name))\
//    var = command_line.next(var);

// some convenient parser abbreviations
#define MAKE_NAME_PARSE_CFG(cfg,type,var,name,def) BASE_PARSE_CFG(cfg,type,var,name,def,static_cast<type>(def))
#define MAKE_PARSE_CFG(cfg,type,var,def) MAKE_NAME_PARSE_CFG(cfg,type,var,#var,def)
#define MAKE_PARSE(type,var,def) MAKE_PARSE_CFG(shellCfgFile,type,var,def)
#define PARSE(var) BASE_PARSE_CFG(shellCfgFile,{},var,#var,var,var)

// from Alessandro. In case the variable already exists
#define PARSE_CLASS_MEMBER(cfg,var,name,def) BASE_PARSE_CFG(cfg,{},var,name,def,def)
#define PARSE_EXTERNAL(cfg,type,var,name,def) \
		extern type var; \
		if (!cfg.have_variable(name) && !command_line.search("-" name))\
		cout << INFO << "WARNING: Parameter " name " is missing! Using default: " << def << std::endl;\
		var = cfg(name, def);\
		if (command_line.search("-" name))\
		var = command_line.next(var);

#define PARSE_EXTERNAL_ARRAY(cfg,type,size,c_var,name,def) \
		extern type* c_var; \
		c_var = new type[size]; \
		for (unsigned int index=0; index<size; ++index) { \
			unsigned int correctIndex=index+1; \
			ostringstream convert; \
			convert<<correctIndex; \
			string stringIndex = convert.str(); \
			string stringName = stringIndex + "/" + name; \
			string optionName = "-" + stringName; \
			if (!cfg.have_variable(stringName) && !command_line.search(optionName)) { \
				cout << INFO << "WARNING: Parameter " << stringName << " is missing! Using default: " << def << std::endl;\
			} \
			c_var[index] = cfg(stringName, def); \
			if (command_line.search(optionName)) { \
				c_var[index] = command_line.next(c_var[index]); \
			} \
		}




// parse an unused variable to prevent it from being added to the UFOs
#define DUMMY_PARSE_CFG(cfg,var) \
		cfg.have_variable(#var);\
		command_line.search("-" #var);

// parse a 3-vector
#define MAKE_VECTOR_PARSE_CFG(cfg,var,def) \
		MAKE_NAME_PARSE_CFG(cfg, std::string, var##_expr, #var, #def "," #def "," #def);\
		mu::Parser p_##var;\
		p_##var.SetExpr(var##_expr);\
		int n_##var;\
		Real * var = p_##var.Eval(n_##var);\
		ASSERT(n_##var == 3)

#define MAKE_VECTOR_PARSE(var,def) MAKE_VECTOR_PARSE_CFG(shellCfgFile,var,def)

// parse a string variable and turn it to uppercase
#define STR_PARSE(var,def) \
		MAKE_NAME_PARSE_CFG(shellCfgFile, std::string, var##_name, #var, #def); \
		std::transform(var##_name.begin(), var##_name.end(), var##_name.begin(), ::toupper);

// read a function from the shell config file and command line, register some variables, and store it in the functions map
#define MAKE_FUNCTION(var,def) \
		MAKE_NAME_PARSE_CFG(shellCfgFile, std::string, var##_expr, #var, #def);\
		Function* var##_func = new Function(var, var##_expr);\
		var##_func->add_variable("t", this->time);\
		var##_func->add_variable("d", density);\
		var##_func->add_variable("x", xyz(0));\
		var##_func->add_variable("y", xyz(1));\
		var##_func->add_variable("z", xyz(2));\
		var##_func->init();\
		functions.insert(std::pair<std::string, Function*>(#var, var##_func));

// function parser for out-of-system evaluation
#define EVAL_BARE_FUNCTION_CFG(cfg,var,def) \
		Real var;\
		MAKE_NAME_PARSE_CFG(cfg, std::string, var##_expr, #var, #def);\
		Function var##_func(var, var##_expr);\
		var##_func.add_variable("t", time);\
		var##_func.add_variable("d", density);\
		var##_func.init();\
		var##_func.evaluate();

#define EVAL_BARE_FUNCTION(var,def) EVAL_BARE_FUNCTION_CFG(shellCfgFile,var,def)

// check if a function is a certain constant c
#define FUNCTION_IS(var,c) (functions[#var]->is_constant() && var == c)

// a shortcut for assigning enums from strings
#define STR_TO_ENUM(var,val) if (var##_name.compare(#val) == 0) var = val;

// measure performance based on compile option (similar to libmesh_logging.h)
#ifdef MEASURE_PERFORMANCE
#include "perf_log.h"
#define PERFLOG_START(name) { perf_log.push(name); }
#define PERFLOG_STOP(name)  { perf_log.pop(name); }
#else
#define PERFLOG_START(name) {}
#define PERFLOG_STOP(name)  {}
#endif

// debugging helpers
#define HERE { cout << INFO << "Here!" << std::endl; }
#define EXIT { cout << INFO << "Exiting!" << std::endl; exit(0); }
#define PRINT(var) { cout << INFO << #var " = " << var << std::endl; }
#define CHECK(var) { if (!std::isfinite(var)) { PRINT(var) exit(1); } }

#endif /* MACROS_H_ */
