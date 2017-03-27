-- Filename: tb_ddfs.vhd
-- Author: Olivier Cotte
-- Date: Jan-2017
-- Description:	assert the relative error above a certain threshold (nominal case is 1e-2)

use work.cordic_types.all;
	
library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use ieee.math_real.all;
use std.textio.all;

library IEEE_PROPOSED;
use IEEE_PROPOSED.FIXED_PKG.ALL;

entity tb_ddfs is
end tb_ddfs;

architecture testbench_pipelined of tb_ddfs is
-- TEST PARAMS --
constant clk_period : time := 1 ns;
constant rst_period : real := 3.5;
constant bins : positive := 10000;
constant num_mode : positive := 2;
constant threshold : real := 0.145;
signal test_cycles_cnt : integer := 0;
signal test_failures : integer := 0;
constant logging : boolean := false;

-- LOG FILE --
constant filename_tb_main : string := "test_results\tb_main.csv";
constant filename_tb_cos : string := "test_results\tb_cos.csv";
constant filename_tb_sin : string := "test_results\tb_sin.csv";
constant filename_tb_arccos : string := "test_results\tb_arccos.csv";
constant filename_tb_arcsin : string := "test_results\tb_arcsin.csv";
constant filename_tb_cosh : string := "test_results\tb_cosh.csv";
constant filename_tb_sinh : string := "test_results\tb_sinh.csv";
constant filename_tb_exp : string := "test_results\tb_exp.csv";
constant filename_tb_sech : string := "test_results\tb_sech.csv";
constant filename_tb_sech_2 : string := "test_results\tb_sech_2.csv";

file res_tb_main : text open write_mode is filename_tb_main;
file res_tb_cos : text open write_mode is filename_tb_cos;
file res_tb_sin : text open write_mode is filename_tb_sin;
file res_tb_arccos : text open write_mode is filename_tb_arccos;
file res_tb_arcsin : text open write_mode is filename_tb_arcsin;
file res_tb_cosh : text open write_mode is filename_tb_cosh;
file res_tb_sinh : text open write_mode is filename_tb_sinh;
file res_tb_exp : text open write_mode is filename_tb_exp;
file res_tb_sech : text open write_mode is filename_tb_sech;
file res_tb_sech_2 : text open write_mode is filename_tb_sech_2;

-- UUT --
signal clk : std_logic;
signal rst_n : std_logic;
signal ddfs_function_type : std_logic_vector(3 downto 0);
signal valid : std_logic;
signal ddfs_function_out : std_logic_vector(BIT_WIDTH-1 downto 0);
signal function_type : unsigned(3 downto 0) := to_unsigned(7, ddfs_function_type'length);
signal cordic_config_i_next, cordic_config_i_ff : cordic_config_t;
signal ddfs_function_i_next, ddfs_function_i_ff : ddfs_function_t;
signal re_init : std_logic := '0';

-- sync with the core
type delay_line_nominal is array (SEND_SYNC_PULSE_PIPELINED_NOMINAL+1 downto 0)of std_logic_vector(BIT_WIDTH-1 downto 0);
type delay_line_hyperbolic_extension_1 is array (SEND_SYNC_PULSE_PIPELINED_EXT_SECH+1 downto 0)of std_logic_vector(BIT_WIDTH-1 downto 0);
type delay_line_hyperbolic_extension_2 is array (SEND_SYNC_PULSE_PIPELINED_EXT_SECH_2+1 downto 0)of std_logic_vector(BIT_WIDTH-1 downto 0);

signal x_nominal_line, y_nominal_line, z_nominal_line : delay_line_nominal := (others => (others => '0'));
signal x_hyperbolic_extension_1_line, y_hyperbolic_extension_1_line, z_hyperbolic_extension_1_line : delay_line_hyperbolic_extension_1 := (others => (others => '0'));
signal x_hyperbolic_extension_2_line, y_hyperbolic_extension_2_line, z_hyperbolic_extension_2_line : delay_line_hyperbolic_extension_2 := (others => (others => '0')); 

signal valid_ff, assert_ff : std_logic := '0';

-- uut input
signal z_in : std_logic_vector(BIT_WIDTH-1 downto 0);

-- trigonometrics functions generator
signal phase_real : real := 0.0;
signal sech_2_real, exp_real, a_cos_h_real, a_sin_h_real : real := 0.0;
signal sech_2_cordic_real, exp_cordic_real, a_cos_h_cordic_real, a_sin_h_cordic_real : real := 0.0;
signal sech_2_diff, exp_diff, a_cos_h_diff, a_sin_h_diff : real := 0.0;
signal sech_2_rel_err, exp_rel_err, a_cos_h_rel_err, a_sin_h_rel_err : real := 0.0;
signal sech_2_min_rel_err, exp_min_rel_err, a_cos_h_min_rel_err, a_sin_h_min_rel_err : real := 1.0;
signal sech_2_max_rel_err, exp_max_rel_err, a_cos_h_max_rel_err, a_sin_h_max_rel_err : real := -1.0;

begin
	------------------------------------------------------------------------------------------------
	-- uut_top_module : 
	------------------------------------------------------------------------------------------------
	uut_top_module : entity work.top_level_ddfs(structural)
	port map(
		clk => clk,
		rst_n => rst_n,
		ddfs_function_type => ddfs_function_type,
		ddfs_function_out => ddfs_function_out
	);
	
	------------------------------------------------------------------------------------------------
	-- uut_state_machine : 
	------------------------------------------------------------------------------------------------
	uut_state_machine : entity work.ddfs_ctrl(rtl)
	port map(
		clk => clk,
		rst_n => rst_n,
		ddfs_function_type => ddfs_function_type,
		cordic_config_sig => cordic_config_i_next,
		cordic_config_reg => cordic_config_i_ff,
		ddfs_function_sig => ddfs_function_i_next,
		ddfs_function_reg => ddfs_function_i_ff,
		valid => valid,
		re_init => re_init
	);
	
	------------------------------------------------------------------------------------------------
	-- uut_phase_gen : 
	------------------------------------------------------------------------------------------------
	uut_phase_gen : entity work.phase_generator(rtl)
	port map (
		clk => clk,
		rst_n => rst_n,
		function_type_sig => ddfs_function_i_next,
		phase_out => z_in
	);
	
	------------------------------------------------------------------------------------------------
	-- tb_function_type_gen : 
	------------------------------------------------------------------------------------------------
	tb_function_type_gen : process
	begin
		ddfs_function_type <= std_logic_vector(function_type);
--		wait;
		function_type <= function_type + 1;
		wait for clk_period * real(bins);
	end process tb_function_type_gen;
	
	------------------------------------------------------------------------------------------------
	-- tb_clk_gen : 
	------------------------------------------------------------------------------------------------
	tb_clk_gen : process
	begin
	    clk <= '1';
	    wait for clk_period / 2;
	    clk <= '0';
	    wait for clk_period / 2;
	end process tb_clk_gen;
	
	------------------------------------------------------------------------------------------------
	-- tb_rst_gen : 
	------------------------------------------------------------------------------------------------
	tb_rst_gen : process
	begin
	    rst_n <= '0';
	    wait for clk_period * rst_period;
	    rst_n <= '1';
		wait;
	    --wait for 100.0 * clk_period * rst_period;
	end process tb_rst_gen;
	
	------------------------------------------------------------------------------------------------
	-- tb_sync_delay_valid : 
	------------------------------------------------------------------------------------------------
	tb_sync_delay_valid : process(clk) is
	begin
		if rising_edge(clk) then
			if rst_n = '0' then
				assert_ff <= '0';
				valid_ff <= '0';
			else
				if valid = '1' then
					valid_ff <= '1';
				elsif re_init = '1' then
					valid_ff <= '0';
				end if;
				assert_ff <= valid_ff; 
			end if;
		end if;
	end process tb_sync_delay_valid;
	
	------------------------------------------------------------------------------------------------
	-- tb_data_in_delay_line : 
	------------------------------------------------------------------------------------------------
	tb_data_in_delay_line : process(clk) is
	begin
		if rising_edge(clk) then
			if rst_n = '0' then
				z_nominal_line <= (others => (others => '0'));
				z_hyperbolic_extension_1_line <= (others => (others => '0'));
				z_hyperbolic_extension_2_line <= (others => (others => '0'));
			else
				case ddfs_function_i_ff is
					when COSINE | SINE | ARCCOSINE | ARCSINE | COSINE_H | SINE_H | EXPONENTIAL =>
					z_nominal_line(z_nominal_line'left-1 downto 0) <= z_nominal_line(z_nominal_line'left downto 1);
					z_nominal_line(z_nominal_line'left) <= z_in;
					
					when SEC_H =>
					z_hyperbolic_extension_1_line(z_hyperbolic_extension_1_line'left-1 downto 0) <= z_hyperbolic_extension_1_line(z_hyperbolic_extension_1_line'left downto 1);
					z_hyperbolic_extension_1_line(z_hyperbolic_extension_1_line'left) <= z_in;
					
					when GEN_SOLITON_SHAPE =>
					z_hyperbolic_extension_2_line(z_hyperbolic_extension_2_line'left-1 downto 0) <= z_hyperbolic_extension_2_line(z_hyperbolic_extension_2_line'left downto 1);
					z_hyperbolic_extension_2_line(z_hyperbolic_extension_2_line'left) <= z_in;
				
					when others =>
				end case;
			end if;
		end if;
	end process tb_data_in_delay_line;
	
	------------------------------------------------------------------------------------------------
	-- tb_validate_output : 
	------------------------------------------------------------------------------------------------
	tb_validate_output : process(clk) is
	variable v_x_in, v_y_in, v_z_in : real;
	variable v_phase_real : real;
	variable v_sech_2_real, v_exp_real, v_a_cos_h_real, v_a_sin_h_real : real;
	variable v_sech_2_cordic_real, v_exp_cordic_real, v_a_cos_h_cordic_real, v_a_sin_h_cordic_real : real;
	variable v_sech_2_diff, v_exp_diff, v_a_cos_h_diff, v_a_sin_h_diff : real;
	variable v_sech_2_rel_err, v_exp_rel_err, v_a_cos_h_rel_err, v_a_sin_h_rel_err : real;
	begin
		if rising_edge(clk) then
			if rst_n = '0' then
				phase_real <= 0.0;
				
				sech_2_real <= 0.0;
				exp_real <= 0.0;
				a_cos_h_real <= 0.0;
				a_sin_h_real <= 0.0;
				
				sech_2_cordic_real <= 0.0;
				exp_cordic_real <= 0.0;
				a_cos_h_cordic_real <= 0.0;
				a_sin_h_cordic_real <= 0.0;
				
				sech_2_diff <= 0.0;
				exp_diff <= 0.0;
				a_cos_h_diff <= 0.0;
				a_sin_h_diff <= 0.0;
				
				sech_2_rel_err <= 0.0;
				exp_rel_err <= 0.0;
				a_cos_h_rel_err <= 0.0;
				a_sin_h_rel_err <= 0.0;
			else
				if valid_ff = '1' then
					case ddfs_function_i_ff is
						when COSINE | ARCCOSINE| COSINE_H  =>
						-- COMPUTE ABSOLUTE AND RELATIVE ERROR
						if ddfs_function_i_ff = COSINE then
							v_z_in := to_real(Circular_Angle_t(z_nominal_line(0)));
							v_a_cos_h_real := cos(v_z_in);	
							v_a_cos_h_cordic_real := to_real(Range_Of_Convergence_t(ddfs_function_out));
						elsif ddfs_function_i_ff = ARCCOSINE then
							v_z_in := to_real(Range_Of_Convergence_t(z_nominal_line(0)));
							v_a_cos_h_real := arccos(v_z_in);	
							v_a_cos_h_cordic_real := to_real(Forbidden_Angle_t(ddfs_function_out));
						else
							v_z_in := to_real(Hyperbolic_Angle_t(z_nominal_line(0)));
							v_a_cos_h_real := cosh(v_z_in);	
							v_a_cos_h_cordic_real := to_real(Hyperbolic_Coordinate_t(ddfs_function_out));
						end if;
						
						v_a_cos_h_diff := (v_a_cos_h_real - v_a_cos_h_cordic_real);
						
						if v_a_cos_h_real = 0.0 then	
							v_a_cos_h_rel_err := abs(v_a_cos_h_diff);
						else			
							v_a_cos_h_rel_err := abs(v_a_cos_h_diff / v_a_cos_h_real);
						end if;
							
						-- UPDATE SIMULATION STATE
						phase_real <= v_z_in;
						a_cos_h_real <= v_a_cos_h_real;
						a_cos_h_cordic_real <= v_a_cos_h_cordic_real;
						a_cos_h_diff <= v_a_cos_h_diff;
						a_cos_h_rel_err <= v_a_cos_h_rel_err;

						when SINE | ARCSINE | SINE_H =>
						-- COMPUTE ABSOLUTE AND RELATIVE ERROR
						if ddfs_function_i_ff = SINE then
							v_z_in := to_real(Circular_Angle_t(z_nominal_line(0)));
							v_a_sin_h_real := sin(v_z_in);	
							v_a_sin_h_cordic_real := to_real(Range_Of_Convergence_t(ddfs_function_out));
						elsif ddfs_function_i_ff = ARCSINE then
							v_z_in := to_real(Range_Of_Convergence_t(z_nominal_line(0)));
							v_a_sin_h_real := arcsin(v_z_in);	
							v_a_sin_h_cordic_real := to_real(Forbidden_Angle_t(ddfs_function_out));
						else
							v_z_in := to_real(Hyperbolic_Angle_t(z_nominal_line(0)));
							v_a_sin_h_real := sinh(v_z_in);	
							v_a_sin_h_cordic_real := to_real(Hyperbolic_Coordinate_t(ddfs_function_out));
						end if;
						
						v_a_sin_h_diff := (v_a_sin_h_real - v_a_sin_h_cordic_real);
						
						if v_a_sin_h_real = 0.0 then	
							v_a_sin_h_rel_err := abs(v_a_sin_h_diff);
						else			
							v_a_sin_h_rel_err := abs(v_a_sin_h_diff / v_a_sin_h_real);
						end if;
							
						-- UPDATE SIMULATION STATE
						phase_real <= v_z_in;
						a_sin_h_real <= v_a_sin_h_real;
						a_sin_h_cordic_real <= v_a_sin_h_cordic_real;
						a_sin_h_diff <= v_a_sin_h_diff;
						a_sin_h_rel_err <= v_a_sin_h_rel_err;

						when EXPONENTIAL =>
						-- COMPUTE ABSOLUTE AND RELATIVE ERROR
						v_z_in := to_real(Hyperbolic_Angle_t(z_nominal_line(0)));
						v_exp_real := exp(v_z_in);	
						v_exp_cordic_real := to_real(Hyperbolic_Coordinate_t(ddfs_function_out));
						v_exp_diff := (v_exp_real - v_exp_cordic_real);
						v_exp_rel_err := abs(v_exp_diff / v_exp_real);
							
						-- UPDATE SIMULATION STATE
						phase_real <= v_z_in;
						exp_real <= v_exp_real;
						exp_cordic_real <= v_exp_cordic_real;
						exp_diff <= v_exp_diff;
						exp_rel_err <= v_exp_rel_err;
						
						when SEC_H | GEN_SOLITON_SHAPE =>
						-- COMPUTE ABSOLUTE AND RELATIVE ERROR
						if ddfs_function_i_ff = SEC_H then
							v_z_in := to_real(Hyperbolic_Angle_t(z_hyperbolic_extension_1_line(0)));
							v_sech_2_real := 1.0 / cosh(v_z_in);	
						else
							v_z_in := to_real(Hyperbolic_Angle_t(z_hyperbolic_extension_2_line(0)));
							v_sech_2_real := 1.0 / cosh(v_z_in);
							v_sech_2_real := v_sech_2_real * v_sech_2_real;
						end if;
						
						-- COMPUTE ABSOLUTE AND RELATIVE ERROR
						v_sech_2_cordic_real := to_real(Range_Of_Convergence_t(ddfs_function_out));
						v_sech_2_diff := (v_sech_2_real - v_sech_2_cordic_real);
						v_sech_2_rel_err := abs(v_sech_2_diff / v_sech_2_real);
							
						-- UPDATE SIMULATION STATE
						phase_real <= v_z_in;
						sech_2_real <= v_sech_2_real;
						sech_2_cordic_real <= v_sech_2_cordic_real;
						sech_2_diff <= v_sech_2_diff;
						sech_2_rel_err <= v_sech_2_rel_err;

						when others =>
					end case;
				end if;
			end if;
		end if;
	end process tb_validate_output;
	
	------------------------------------------------------------------------------------------------
	-- tb_assert_output : 
	------------------------------------------------------------------------------------------------
	tb_assert_output : process(clk) is
	begin
		if rising_edge(clk) then
			if rst_n = '0' then
			else
				if assert_ff = '1' then
					case ddfs_function_i_ff is
						when COSINE | COSINE_H | ARCCOSINE =>
						assert a_cos_h_rel_err < threshold
						report real'image(a_cos_h_rel_err) & " > " & real'image(threshold) severity error;
						if a_cos_h_rel_err > threshold then
							test_failures <= test_failures + 1;
						end if;
						if a_cos_h_rel_err < a_cos_h_min_rel_err then
							a_cos_h_min_rel_err <= a_cos_h_rel_err;
						end if;
						if a_cos_h_rel_err > a_cos_h_max_rel_err then
							a_cos_h_max_rel_err <= a_cos_h_rel_err;
						end if;

						when SINE | SINE_H | ARCSINE =>
						assert a_sin_h_rel_err < threshold 
						report real'image(a_sin_h_rel_err) & " > " & real'image(threshold) severity error;
						if a_sin_h_rel_err > threshold then
							test_failures <= test_failures + 1;
						end if;
						if a_sin_h_rel_err < a_sin_h_min_rel_err then
							a_sin_h_min_rel_err <= a_sin_h_rel_err;
						end if;
						if a_sin_h_rel_err > a_sin_h_max_rel_err then
							a_sin_h_max_rel_err <= a_sin_h_rel_err;
						end if;
						
						when EXPONENTIAL =>
						assert exp_rel_err < threshold 
						report real'image(exp_rel_err) & " > " & real'image(threshold) severity error;
						if exp_rel_err > threshold then
							test_failures <= test_failures + 1;
						end if;
						if exp_rel_err < exp_min_rel_err then
							exp_min_rel_err <= exp_rel_err;
						end if;
						if exp_rel_err > exp_max_rel_err then
							exp_max_rel_err <= exp_rel_err;
						end if;
						
						when SEC_H | GEN_SOLITON_SHAPE =>
						assert sech_2_rel_err < threshold 
						report real'image(sech_2_rel_err) & " > " & real'image(threshold) severity error;
						if sech_2_rel_err > threshold then
							test_failures <= test_failures + 1;
						end if;
						if sech_2_rel_err < sech_2_min_rel_err then
							sech_2_min_rel_err <= sech_2_rel_err;
						end if;
						if sech_2_rel_err > sech_2_max_rel_err then
							sech_2_max_rel_err <= sech_2_rel_err;
						end if;
						
						when others =>
					end case;
				end if;
			end if;
		end if;
	end process tb_assert_output;
	
	------------------------------------------------------------------------------------------------
	-- tb_cycles_count : 
	------------------------------------------------------------------------------------------------
	tb_cycles_count : process(clk) is
	variable out_line : line;
	begin
		if rising_edge(clk) then
			if test_cycles_cnt < (num_mode * bins) then
				test_cycles_cnt <= test_cycles_cnt + 1;
			else
				write(out_line, "test_failures" & "," & integer'image(test_failures));
				writeline(res_tb_main, out_line);
				assert false report "endsim" severity failure;
 				--WAIT; --allows the simulation to halt!
			end if;
		end if;
	end process tb_cycles_count;
	
	------------------------------------------------------------------------------------------------
	-- init_file : 
	------------------------------------------------------------------------------------------------
	init_file : process
	variable out_line_tb_cos : line;
	variable out_line_tb_sin : line;
	variable out_line_tb_arccos : line;
	variable out_line_tb_arcsin : line;
	variable out_line_tb_cosh : line;
	variable out_line_tb_sinh : line;
	variable out_line_tb_exp : line;
	variable out_line_tb_sech : line;
	variable out_line_tb_sech_2 : line;
	begin
		if logging = true then
			write(out_line_tb_cos,
			"phase(rad)" &
			",value_cordic" &
			",real_value" &
			",abs_err" &
			",rel_err");
			
		    write(out_line_tb_cos,
			"phase(rad)" &
			",cos_cordic" &
			",cos_real" &
			",abs_err_cos" &
			",rel_err_cos");
					
			write(out_line_tb_sin,
			"phase(rad)" &
			",sin_cordic" &
			",sin_real" &
			",abs_err_sin" &
			",rel_err_sin");
					
			write(out_line_tb_arccos,
			"phase(rad)" &
			",acos_cordic" &
			",acos_real" &
			",abs_err_acos" &
			",rel_err_acos");
				
			write(out_line_tb_arcsin,
			"phase(rad)" &
			",asin_cordic" &
			",asin_real" &
			",abs_err_asin" &
			",rel_err_asin");
					
			write(out_line_tb_cosh,
			"phase(rad)" &
			",cosh_cordic" &
			",cosh_real" &
			",abs_err_cosh" &
			",rel_err_cosh");
					
			write(out_line_tb_sinh,
			"phase(rad)" &
			",sinh_cordic" &
			",sinh_real" &
			",abs_err_sinh" &
			",rel_err_sinh");
					
			write(out_line_tb_exp,
			"phase(rad)" &
			",exp_cordic" &
			",exp_real" &
			",abs_err_exp" &
			",rel_err_exp");
					
			write(out_line_tb_sech,
			"phase(rad)" &
			",sech_cordic" &
			",sech_real" &
			",abs_err_sech" &
			",rel_err_sech");
				
			write(out_line_tb_sech_2,
			"phase(rad)" &
			",sech_2_cordic" &
			",sech_2_real" &
			",abs_err_sech_2" &
			",rel_err_sech_2");
			
			writeline(res_tb_cos, out_line_tb_cos);
			writeline(res_tb_sin, out_line_tb_sin);
			writeline(res_tb_arccos, out_line_tb_arccos);
			writeline(res_tb_arcsin, out_line_tb_arcsin);
			writeline(res_tb_cosh, out_line_tb_cosh);
			writeline(res_tb_sinh, out_line_tb_sinh);
			writeline(res_tb_exp, out_line_tb_exp);
			writeline(res_tb_sech, out_line_tb_sech);
			writeline(res_tb_sech_2, out_line_tb_sech_2);
		end if;
	    wait;
	end process init_file;
	
	------------------------------------------------------------------------------------------------
	-- log_to_file : 
	------------------------------------------------------------------------------------------------
	log_to_file : process(clk)
	variable out_line : line;
	begin
		if rising_edge(clk) then
			if logging = true and assert_ff = '1' then
				case ddfs_function_i_ff is
					when COSINE | ARCCOSINE | COSINE_H =>
					if a_cos_h_rel_err < threshold then
						write(out_line, real'image(phase_real) &
						"," & real'image(a_cos_h_real) &
						"," & real'image(a_cos_h_cordic_real) &
						"," & real'image(a_cos_h_diff) &
						"," & real'image(a_cos_h_rel_err));
					else
						write(out_line, real'image(phase_real) &
						"," & real'image(a_cos_h_real) &
						"," & real'image(a_cos_h_cordic_real) &
						"," & real'image(a_cos_h_diff) &
						"," & real'image(a_cos_h_rel_err) & ",*");
					end if;
					
					writeline(res_tb_main, out_line);
					if ddfs_function_i_ff = COSINE then
						writeline(res_tb_cos, out_line);
					elsif ddfs_function_i_ff = ARCCOSINE then
						writeline(res_tb_arccos, out_line);
					else
						writeline(res_tb_cosh, out_line);
					end if;
					
					when SINE | ARCSINE | SINE_H =>
					if a_sin_h_rel_err < threshold then
						write(out_line, real'image(phase_real) &
						"," & real'image(a_sin_h_real) &
						"," & real'image(a_sin_h_cordic_real) &
						"," & real'image(a_sin_h_diff) &
						"," & real'image(a_sin_h_rel_err));
					else
						write(out_line, real'image(phase_real) &
						"," & real'image(a_sin_h_real) &
						"," & real'image(a_sin_h_cordic_real) &
						"," & real'image(a_sin_h_diff) &
						"," & real'image(a_sin_h_rel_err) & ",*");
					end if;
					
					writeline(res_tb_main, out_line);
					if ddfs_function_i_ff = SINE then
						writeline(res_tb_sin, out_line);
					elsif ddfs_function_i_ff = ARCSINE then
						writeline(res_tb_arcsin, out_line);
					else
						writeline(res_tb_sinh, out_line);
					end if;
					
					when EXPONENTIAL =>
					if exp_rel_err < threshold then
						write(out_line, real'image(phase_real) &
						"," & real'image(exp_real) &
						"," & real'image(exp_cordic_real) &
						"," & real'image(exp_diff) &
						"," & real'image(exp_rel_err));
					else
						write(out_line, real'image(phase_real) &
						"," & real'image(exp_real) &
						"," & real'image(exp_cordic_real) &
						"," & real'image(exp_diff) &
						"," & real'image(exp_rel_err) & ",*");
					end if;
					writeline(res_tb_main, out_line);
					writeline(res_tb_exp, out_line);
						
					when SEC_H | GEN_SOLITON_SHAPE =>
						if sech_2_rel_err < threshold then
						write(out_line, real'image(phase_real) &
						"," & real'image(sech_2_real) &
						"," & real'image(sech_2_cordic_real) &
						"," & real'image(sech_2_diff) &
						"," & real'image(sech_2_rel_err));
					else
						write(out_line, real'image(phase_real) &
						"," & real'image(sech_2_real) &
						"," & real'image(sech_2_cordic_real) &
						"," & real'image(sech_2_diff) &
						"," & real'image(sech_2_rel_err) & ",*");
					end if;
					
					writeline(res_tb_main, out_line);
					if ddfs_function_i_ff = SEC_H then
						writeline(res_tb_sech, out_line);
					else
						writeline(res_tb_sech_2, out_line);
					end if;

					when others =>
				end case;
			end if;
		end if;
	end process log_to_file;
end testbench_pipelined;