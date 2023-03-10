// $Id: Solvers.cpp,v 1.49.2.1 2006-03-10 15:54:18 geuzaine Exp $
//
// Copyright (C) 1997-2006 C. Geuzaine, J.-F. Remacle
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA.
// 
// Please report all bugs and problems to <gmsh@geuz.org>.

#include "Gmsh.h"
#include "Solvers.h"
#include "GmshServer.h"
#include "OpenFile.h"
#include "GmshUI.h"
#include "GUI.h"
#include "Mesh.h"
#include "Draw.h"
#include "Context.h"

extern Context_T CTX;
extern GUI *WID;

SolverInfo SINFO[MAXSOLVERS];

int GmshServer::init = 0;
int GmshServer::s;

// This routine polls the socket at least every 'waitint' seconds and
// returns 0 if data is available or 1 if there was en error or if the
// process was killed. Otherwise it just tends to current GUI events
// (this is easier to manage than non-blocking IO, and simpler than
// using the "real" solution, i.e., threads. Another possibility would
// be to use Fl::add_fd())

int WaitForData(int socket, int num, double waitint)
{
  while(1){
    if((num >= 0 && SINFO[num].pid < 0) || (num < 0 && !CTX.solver.listen)){
      // process has been killed or we stopped listening
      return 1;
    }

    // check if there is data (call select with a zero timeout to
    // return immediately, i.e., do polling)
    int ret = myselect(socket, 0);

    if(ret == 0){ 
      // nothing available: wait at most waitint seconds
      WID->wait(waitint);
    }
    else if(ret > 0){ 
      // data is there
      return 0;
    }
    else{ 
      // an error happened
      if(num >= 0)
	SINFO[num].pid = -1;
      return 1;
    }
  }
}

// This routine either launches a solver and waits for some answer (if
// num >= 0), or simply waits for messages (if num < 0)

int Solver(int num, char *args)
{
  char command[1024], sockname[1024], prog[1024], tmp[1024], tmp2[1024];

 new_connection:

  GmshServer server(CTX.solver.max_delay);

  if(num >= 0){
    FixWindowsPath(SINFO[num].executable_name, prog);
    if(!SINFO[num].client_server) {
      sprintf(command, "%s %s", prog, args);
#if !defined(WIN32)
      strcat(command, " &");
#endif
      server.StartClient(command);
      return 1;
    }
  }
  else{
    if(!CTX.solver.listen){
      Msg(INFO, "Stopped listening for solver connections");
      return 0;
    }
    // we don't know who will (maybe) contact us
    strcpy(prog, "");
    strcpy(command, "");
  }

  if(!strstr(CTX.solver.socket_name, ":")){
    // Unix socket
    if(num >= 0)
      sprintf(tmp, "%s%s-%d", CTX.home_dir, CTX.solver.socket_name, num);
    else
      sprintf(tmp, "%s%s", CTX.home_dir, CTX.solver.socket_name);
    FixWindowsPath(tmp, sockname);
  }
  else{
    // TCP/IP socket
    strcpy(sockname, CTX.solver.socket_name);
  }

  if(num >= 0){
    sprintf(tmp, "\"%s\"", sockname);
    sprintf(tmp2, SINFO[num].socket_command, tmp);
    sprintf(command, "%s %s %s", prog, args, tmp2);
#if !defined(WIN32)
    strcat(command, " &");
#endif
  }

  int sock = server.StartClient(command, sockname);

  if(sock < 0) {
    switch (sock) {
    case -1:
      Msg(GERROR, "Couldn't create socket '%s'", sockname);
      break;
    case -2:
      Msg(GERROR, "Couldn't bind socket to name '%s'", sockname);
      break;
    case -3:
      Msg(GERROR, "Socket listen failed on '%s'", sockname);
      break;
    case -4:
      Msg(GERROR, "Socket listen timeout on '%s'", sockname);
      Msg(GERROR, "Is '%s' correctly installed?", prog);
      break;
    case -5:
      Msg(GERROR, "Socket accept failed on '%s'", sockname);
      break;
    case -6:
      Msg(INFO, "Stopped listening for solver connections");
      server.StopClient();
      break;
    case -7:
      Msg(GERROR, "Unix sockets not available on Windows without Cygwin");
      Msg(GERROR, "Use TCP/IP sockets instead");
      break;
    case -8:
      Msg(GERROR, "Could not initialize Windows sockets");
      break;
    }
    if(num >= 0){
      for(int i = 0; i < SINFO[num].nboptions; i++)
	WID->solver[num].choice[i]->clear();
    }
    return 0;
  }

  if(num >= 0){
    for(int i = 0; i < SINFO[num].nboptions; i++)
      SINFO[num].nbval[i] = 0;
    SINFO[num].pid = 0;
  }

  while(1) {

    int stop = (num < 0 && !CTX.solver.listen);

    if(stop || (num >= 0 && SINFO[num].pid < 0))
      break;

    stop = WaitForData(sock, num, 0.1);

    if(stop || (num >= 0 && SINFO[num].pid < 0))
      break;

    int type, length;
    if(server.ReceiveMessageHeader(&type, &length)){
      char *message = new char[length + 1];
      if(server.ReceiveMessageBody(length, message)){
	switch (type) {
	case GmshServer::CLIENT_START:
	  if(num >= 0)
	    SINFO[num].pid = atoi(message);
	  break;
	case GmshServer::CLIENT_STOP:
	  stop = 1;
	  if(num >= 0)
	    SINFO[num].pid = -1;
	  break;
	case GmshServer::CLIENT_PROGRESS:
	  if(num >= 0)
	    Msg(STATUS3N, "%s %s", SINFO[num].name, message);
	  else
	    Msg(STATUS3N, "%s", message);
	  break;
	case GmshServer::CLIENT_OPTION_1:
	  if(num >= 0)
	    strcpy(SINFO[num].option[0][SINFO[num].nbval[0]++], message);
	  break;
	case GmshServer::CLIENT_OPTION_2:
	  if(num >= 0)
	    strcpy(SINFO[num].option[1][SINFO[num].nbval[1]++], message);
	  break;
	case GmshServer::CLIENT_OPTION_3:
	  if(num >= 0)
	    strcpy(SINFO[num].option[2][SINFO[num].nbval[2]++], message);
	  break;
	case GmshServer::CLIENT_OPTION_4:
	  if(num >= 0)
	    strcpy(SINFO[num].option[3][SINFO[num].nbval[3]++], message);
	  break;
	case GmshServer::CLIENT_OPTION_5:
	  if(num >= 0)
	    strcpy(SINFO[num].option[4][SINFO[num].nbval[4]++], message);
	  break;
	case GmshServer::CLIENT_MERGE_FILE:
	  if(num < 0 || (num >= 0 && SINFO[num].merge_views)) {
	    int n = List_Nbr(CTX.post.list);
	    MergeProblem(message);
	    Draw();
	    if(n != List_Nbr(CTX.post.list))
	      WID->set_context(menu_post, 0);
	  }
	  break;
	case GmshServer::CLIENT_PARSE_STRING:
	  ParseString(message);
	  Draw();
	  break;
	case GmshServer::CLIENT_INFO:
	  Msg(SOLVER, "%-8.8s: %s", num >= 0 ? SINFO[num].name : "Client", message);
	  break;
	case GmshServer::CLIENT_WARNING:
	case GmshServer::CLIENT_ERROR:
	  Msg(SOLVERR, "%-8.8s: %s", num >= 0 ? SINFO[num].name : "Client", message);
	  break;
	default:
	  Msg(WARNING, "Unknown type of message received from %s",
	      num >= 0 ? SINFO[num].name : "client");
	  Msg(SOLVER, "%-8.8s: %s", num >= 0 ? SINFO[num].name : "Client", message);
	  break;
	}
	WID->check();
      }
      else{
	Msg(WARNING, "Failed to receive message body on socket: aborting");
	break;
      }
      delete [] message;
    }
    else{
      // didn't get any header, just abort
      break;
    }
  }
  
  if(num >= 0){
    for(int i = 0; i < SINFO[num].nboptions; i++) {
      if(SINFO[num].nbval[i]) {
	WID->solver[num].choice[i]->clear();
	for(int j = 0; j < SINFO[num].nbval[i]; j++)
	  WID->solver[num].choice[i]->add(SINFO[num].option[i][j]);
	WID->solver[num].choice[i]->value(0);
      }
    }
  }

  if(server.StopClient() < 0)
    Msg(WARNING, "Impossible to unlink the socket '%s'", sockname);

  if(num >= 0){
    Msg(STATUS3N, "Ready");
  }
  else{
    Msg(INFO, "Client disconnected: starting new connection");
    goto new_connection;
  }

  return 1;
}
