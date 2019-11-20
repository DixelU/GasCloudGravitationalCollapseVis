#pragma once
#define SMP_BOOL_SETTINGS_EMPTY_TRACKS_RMV 0b1u
#define SMP_BOOL_SETTINGS_REMNANTS_RMV 0b10u
#define SMP_BOOL_SETTINGS_ALL_INSTRUMENTS_TO_PIANO 0b100u
#define SMP_BOOL_SETTINGS_IGNORE_TEMPOS 0b1000u
#define SMP_BOOL_SETTINGS_IGNORE_ALL_BUT_TEMPOS_NOTES_AND_PITCH 0b10000u
#define SMP_BOOL_SETTINGS_IGNORE_PITCHES 0b100000u
#define SMP_BOOL_SETTINGS_IGNORE_NOTES 0b1000000u

enum _BoolSettings {
	remove_empty_tracks = SMP_BOOL_SETTINGS_EMPTY_TRACKS_RMV,
	remove_remnants = SMP_BOOL_SETTINGS_REMNANTS_RMV,
	all_instruments_to_piano = SMP_BOOL_SETTINGS_ALL_INSTRUMENTS_TO_PIANO,
	ignore_tempos = SMP_BOOL_SETTINGS_IGNORE_TEMPOS,
	ignore_all_but_tempos_notes_and_pitch = SMP_BOOL_SETTINGS_IGNORE_ALL_BUT_TEMPOS_NOTES_AND_PITCH,
	ignore_pitches = SMP_BOOL_SETTINGS_IGNORE_PITCHES,
	ignore_notes = SMP_BOOL_SETTINGS_IGNORE_NOTES
};

#define TT_UNSPECIFIED 0b0
#define TT_INPUT_FIELD 0b1
#define TT_MOVEABLE_WINDOW 0b10
#define TT_BUTTON 0b100
#define TT_TEXTBOX 0b1000
#define TT_SELPROPLIST 0b10000
#define TT_CHECKBOX 0b100000

enum _TellType {
	unspecified = TT_UNSPECIFIED,
	input_field = TT_INPUT_FIELD,
	moveable_window = TT_MOVEABLE_WINDOW,
	button = TT_BUTTON,
	textbox = TT_TEXTBOX,
	selectable_properted_list = TT_SELPROPLIST,
	checkbox = TT_CHECKBOX
};

#define GLOBAL_LEFT 0b0001
#define GLOBAL_RIGHT 0b0010
#define GLOBAL_TOP 0b0100
#define GLOBAL_BOTTOM 0b1000

enum _Positioning {
	vertical=0b1,
	horizontal=0b10
};
enum _Align{
	center=0,
	left=GLOBAL_LEFT,
	right=GLOBAL_RIGHT,
	top=GLOBAL_TOP,
	bottom=GLOBAL_BOTTOM
};
const double pi = atan(1) * 4;