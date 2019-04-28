/*************************************************************************
	> File Name: keyboard.h
	> Author: 
	> Mail: 
	> Created Time: 2019年04月16日 星期二 11时53分24秒
 ************************************************************************/

#ifndef _KEYBOARD_H
#define _KEYBOARD_H
#endif

#include <stdio.h>
#include <string.h>
#include <termios.h>

void init_keyboard(void);
void close_keyboard(void);
int kbhit(void);
int readch(void);

