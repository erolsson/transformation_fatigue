        SUBROUTINE PART_TYPE(LOCNUM, SIDE)
        CHARACTER NAME*255
        EXTERNAL GETPARTINFO
        CALL GETPARTINFO(LOCNUM, JTYPE, NAME, USER_NUMBER, JRCD)
        END

        SUBROUTINE GETPARTINFOC(NAME, LOCNUM, JTYPE, USER_NUMBER, JRCD)
        EXTERNAL GETPARTINFO
        CHARACTER NAME*80
        CALL GETPARTINFO(LOCNUM, JTYPE, NAME, USER_NUMBER, JRCD)
        END

        SUBROUTINE GETELEMNUMBERUSER(LOCNUM, USER_NUMBER)
        EXTERNAL GETPARTINFO
        CHARACTER NAME*80
        CALL GETPARTINFO(LOCNUM, 1, NAME, USER_NUMBER, JRCD)
        END

        SUBROUTINE GETPARTNAME(LOCNUM, NAME)
        EXTERNAL GETPARTINFO
        CHARACTER NAME*80
        CALL GETPARTINFO(LOCNUM, 1, NAME, USER_NUMBER, JRCD)
        END