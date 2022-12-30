typedef union {
  char *c;
  int i;
  unsigned int u;
  double d;
  double v[5];
  Shape s;
  List_T *l;
} YYSTYPE;
#define	tDOUBLE	257
#define	tSTRING	258
#define	tBIGSTR	259
#define	tEND	260
#define	tAFFECT	261
#define	tDOTS	262
#define	tPi	263
#define	tMPI_Rank	264
#define	tMPI_Size	265
#define	tExp	266
#define	tLog	267
#define	tLog10	268
#define	tSqrt	269
#define	tSin	270
#define	tAsin	271
#define	tCos	272
#define	tAcos	273
#define	tTan	274
#define	tRand	275
#define	tAtan	276
#define	tAtan2	277
#define	tSinh	278
#define	tCosh	279
#define	tTanh	280
#define	tFabs	281
#define	tFloor	282
#define	tCeil	283
#define	tFmod	284
#define	tModulo	285
#define	tHypot	286
#define	tPrintf	287
#define	tSprintf	288
#define	tStrCat	289
#define	tStrPrefix	290
#define	tStrRelative	291
#define	tBoundingBox	292
#define	tDraw	293
#define	tToday	294
#define	tPoint	295
#define	tCircle	296
#define	tEllipse	297
#define	tLine	298
#define	tSurface	299
#define	tSpline	300
#define	tVolume	301
#define	tCharacteristic	302
#define	tLength	303
#define	tParametric	304
#define	tElliptic	305
#define	tPlane	306
#define	tRuled	307
#define	tTransfinite	308
#define	tComplex	309
#define	tPhysical	310
#define	tUsing	311
#define	tBump	312
#define	tProgression	313
#define	tPlugin	314
#define	tRotate	315
#define	tTranslate	316
#define	tSymmetry	317
#define	tDilate	318
#define	tExtrude	319
#define	tDuplicata	320
#define	tLoop	321
#define	tRecombine	322
#define	tDelete	323
#define	tCoherence	324
#define	tIntersect	325
#define	tAttractor	326
#define	tLayers	327
#define	tAlias	328
#define	tAliasWithOptions	329
#define	tText2D	330
#define	tText3D	331
#define	tInterpolationScheme	332
#define	tTime	333
#define	tCombine	334
#define	tBSpline	335
#define	tBezier	336
#define	tNurbs	337
#define	tOrder	338
#define	tWith	339
#define	tBounds	340
#define	tKnots	341
#define	tColor	342
#define	tColorTable	343
#define	tFor	344
#define	tIn	345
#define	tEndFor	346
#define	tIf	347
#define	tEndIf	348
#define	tExit	349
#define	tReturn	350
#define	tCall	351
#define	tFunction	352
#define	tTrimmed	353
#define	tShow	354
#define	tHide	355
#define	tGetValue	356
#define	tGMSH_MAJOR_VERSION	357
#define	tGMSH_MINOR_VERSION	358
#define	tGMSH_PATCH_VERSION	359
#define	tAFFECTPLUS	360
#define	tAFFECTMINUS	361
#define	tAFFECTTIMES	362
#define	tAFFECTDIVIDE	363
#define	tOR	364
#define	tAND	365
#define	tEQUAL	366
#define	tNOTEQUAL	367
#define	tAPPROXEQUAL	368
#define	tLESSOREQUAL	369
#define	tGREATEROREQUAL	370
#define	tCROSSPRODUCT	371
#define	tPLUSPLUS	372
#define	tMINUSMINUS	373
#define	UNARYPREC	374


extern YYSTYPE yylval;