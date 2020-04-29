Module Graph1

    'Construction de graphique
    Public Sub manage_graph(ByRef frm02 As Chloride_transport.frmChlor, ByRef frm03 As Chloride_transport.frmScale1)
        Dim cpteur As Short = 0
        Dim Filtre As String
        Dim Index As Short
        Dim Directoire As Boolean
        Dim Titre As String
        Dim Canc As Boolean = False
        Dim OutFile As String
        Dim a As Single
        Dim NoAxeX As Short
        Dim NoAxeY As Short
        Dim EctY As Single
        Dim EctX As Single
        Dim sh, ih, ix, sx As Single
        Dim hEcart As Single
        Dim xEcart As Single

        a = frm02.Width / 34

        Dim RFile As String
        Dim TestB As String
        Dim iPos As Short
        Dim iPos_old As Short
        Dim Boucle1 As Short
        Dim Boucle2 As Short
        Dim UserInput2 As String
        Dim PostFile As String
        Dim NCurve As Short
        Dim NTCurve As Short
        Dim nTProfil As Short
        Dim nrofil As Short
        Dim TesGraph As Short
        Dim Data(1, 1) As Object
        Dim XMI, XMA, YMI, YMA As Single
        Dim Titel, LibelX, LibelY As String
        Dim x1, y1, x2, y2 As Short
        Dim CoorX As Single
        Dim Dim1, Nline, kk, j As Short
        Dim Color As Color
        Dim i As Short
        Dim TesPoint As Short = 0
        Dim x1Min As Single
        Dim x1Max As Single
        Dim x2Min As Single
        Dim x2Max As Single
        Dim yMin As Single
        Dim yMax As Single
        Dim NCol As Short
        Dim UserInput1 As Short

        frm02.Command1.Visible = False
        frm02.PictureBox6.Visible = False
        'R�cup�rer un nombre de l'utilisateur

        UserInput2 = InputBox("Entrez le nombre total de courbes :", "graphique n�" & Boucle1, CStr(1))
        Msg_noEntry(UserInput2, Canc)
        If Canc = False Then Msg_noNumeric(UserInput2, Canc)
        If Canc = True Then Exit Sub
        NTCurve = CShort(UserInput2)

        TesGraph = 0
        For Boucle2 = 1 To NTCurve
            NCurve = Boucle2
            If TesGraph = 1 Then
                TesGraph = MsgBox("Utilisation de l'ancien fichier de donn�es ?", MsgBoxStyle.YesNo, "Fichier de donn�es")
                If TesGraph = MsgBoxResult.Yes Then
                    TesGraph = 1
                Else
                    TesGraph = 2
                End If
            End If
            read_file_result(RFile, Data, TesGraph, Nline, NCol, NCurve, x2Min, x2Max, YMA, YMI, x1Min, x1Max, yMin, yMax, Canc)
            If Canc = True Then Exit Sub
            Design_graph(frm02, cpteur, TesGraph, Data, NCurve, NoAxeX, NoAxeY, RFile, XMI, XMA, YMI, YMA, NCol, Nline, x2Min, x2Max, x1Min, x1Max, yMin, yMax)
        Next Boucle2
        Kill("tttemp.txt")
        'R�cup�rer un nombre de l'utilisateur
        Try
            UserInput1 = InputBox("Entrez le nombre total de profils mesur�s :", "Editeur de graphique", CStr(0))
        Catch ex As Exception
            GoTo c
        End Try
        'Msg_noEntry(UserInput1, Canc)
        'If Canc = False Then Msg_noNumeric(UserInput1, Canc)
        'If Canc = True Then Exit Sub

        nTProfil = CShort(UserInput1)
        For Boucle2 = 1 To nTProfil
            nrofil = Boucle2
            Design_point(frm02, Data, nrofil, NoAxeX, NoAxeY, TesPoint, XMI, XMA, YMI, YMA)
        Next Boucle2

c:      Scale_info(frm02, frm03, Data, NTCurve, nTProfil, TesGraph, NoAxeX, NoAxeY, a, TesPoint, XMI, XMA, YMI, YMA, EctX, EctY)

        GraphBMP(frm02, nTProfil, NTCurve, Nline, CoorX, kk, NoAxeX, NoAxeY, XMI, XMA, YMI, YMA, EctX, EctY)
        TestB = MsgBox("Changer l'�chelle du graphique ?", MsgBoxStyle.YesNo, "Fichier de donn�es")
        If TestB = MsgBoxResult.Yes Then GoTo c

        ''''''''''''''''''''''''''''''''''
        Filtre = "txt files (GR_*.BMP)|GR_*.BMP"
        Index = 1
        Directoire = True
        Titre = "Enregistrer le graphique"
        SaveDialog(OutFile, Canc, Filtre, Index, Directoire, Titre)
        If Canc = True Then GoTo b

        ''''''''''''''''''''''''''''''''''
        FilePost(OutFile, PostFile)
        FileOnly(OutFile)
        OutFile = PostFile & "GR_" & OutFile
        frm02.PictureBox1.Image.Save(OutFile)

b:
        For Boucle2 = 1 To NTCurve
            OutFile = "G_temp_" & Boucle2 & ".txt"
            Kill(OutFile)
        Next Boucle2
        For Boucle2 = 1 To nTProfil
            OutFile = "P_temp_" & Boucle2 & ".txt"
            Kill(OutFile)
        Next Boucle2
        frm02.PictureBox6.Visible = True
        frm02.Hide()
        frm02.Close()

    End Sub

    'Lecture du fichier de r�sultats
    Public Sub read_file_result(ByRef RFile As String, ByRef Data As Object, ByRef TesGraph As Short, ByRef Nline As Short, ByRef Ncol As Short, ByRef NCurve As Short, ByRef x2Min As Single, ByRef x2Max As Single, ByRef YMA As Single, ByRef YMI As Single, ByRef x1Min As Single, ByRef x1Max As Single, ByRef yMin As Single, ByRef yMax As Single, ByRef Canc As Boolean)
        Dim k As Short
        Dim Dim2 As String
        Dim iPos As Short
        Dim nFic As Short
        Dim Dim1 As String
        Dim i As Short
        Dim j As Short
        Dim Filtre As String
        Dim Index As Short
        Dim Directoire As Boolean
        Dim Titre As String
        Canc = False
        Dim OutFile As String

        If TesGraph >= 0 Or TesGraph = 2 Then ' test si on active la fonction depuis compute_all ou depuis mdi_chlor
            If TesGraph <> 1 Then
                ''''''''''''''''''''''''''''''''''''''''''''
                Filtre = "txt files (R_*.txt)|R_*.txt"
                Index = 2
                Directoire = True
                Titre = "Ouvrir le fichier r�sultat, courbe n� " & NCurve
                OpenDialog(RFile, Canc, Filtre, Index, Directoire, Titre)
                If Canc = True Then GoTo b
                '''''''''''''''''''''''''''''''''''''''''''''
            End If
            nFic = FreeFile()
            FileOpen(nFic, RFile, OpenMode.Input, OpenAccess.Read, OpenShare.Shared)
            OutFile = RFile
            FileOnly(OutFile)

            Input(nFic, Dim1) 'nbre de ligne
            Input(nFic, Dim2) 'nbre de colonnes
            Dim1 = "_"
            iPos = InStr(1, Dim2, Dim1, CompareMethod.Text) 'position de "_"
            Nline = CShort(Mid(Dim2, 1, iPos - 1))
            Dim2 = CStr(Mid(Dim2, iPos + 1))
            Dim1 = " "
            iPos = 1
            Do While iPos = 1
                iPos = InStr(1, Dim2, Dim1, CompareMethod.Text)
                Dim2 = CStr(Mid(Dim2, iPos + 1))
                iPos = InStr(1, Dim2, Dim1, CompareMethod.Text)
            Loop
            Ncol = CShort(Mid(Dim2, 1, iPos + 1))
            ReDim Data(Nline + 2, Ncol)
            For i = 3 To Ncol
                Input(nFic, Data(1, i))
            Next i

            For j = 2 To Nline + 2
                For k = 1 To Ncol
                    Try
                        Input(nFic, Data(j, k)) 'lecture du fichier r�sultat
                    Catch ex As Exception When Data(j, k) = ""
                        Exit For
                        Exit For
                    End Try
                Next k
            Next j
            FileClose(nFic)

            'valeur des axes max et min
            x1Min = Data(1, 3) 'profondeur
            If CType(Data(1, Ncol - 1), String) = "cumul" Then
                x1Max = Data(1, Ncol - 2)
            Else
                x1Max = Data(1, Ncol - 1)
            End If
            x2Min = Data(2, 2) 'temps
            x2Max = Data(Nline + 2, 2)

            yMin = 9999999      'r�sultats
            yMax = 0
            For j = 2 To Nline + 2
                For k = 3 To Ncol - 1
                    If Data(j, k) < yMin Then yMin = Data(j, k)
                    If Data(j, k) > yMax Then yMax = Data(j, k)
                Next k
            Next j
            If TesGraph = 0 Then
                YMA = yMax
                YMI = yMin
            End If

        Else
            'si depuis compute_all
        End If

b:      'user pressed cancel error

    End Sub

    'Cherche le fichier r�sultat, g�re le choix des axes, pointe la ligne ou la colonne de r�sultat et enregistre dans un fichier texte
    Public Sub Design_graph(ByRef frm02 As Chloride_transport.frmChlor, ByRef cpteur As Short, ByRef TesGraph As Short, ByRef Data As Object, ByRef NCurve As Short, ByRef NoAxeX As Short, ByRef NoAxeY As Short, ByRef RFile As String, ByRef XMI As Single, ByRef XMA As Single, ByRef YMI As Single, ByRef YMA As Single, ByRef NCol As Short, ByRef NLine As Short, ByRef x2Min As Single, ByRef x2Max As Single, ByRef x1Min As Single, ByRef x1Max As Single, ByRef yMin As Single, ByRef yMax As Single)
        Dim CoorX As Single
        Dim TypY As String
        Dim i As Short
        Dim k As Short
        Dim CoPosX As Short
        Dim nFic As Short
        Dim outfile As String
        Dim Rfileprov As String

        Dim frm01 As New frmAxe

        If TesGraph = 0 Then
            'choix pour l'axe des x entre la profondeur et le temps
            Rfileprov = RFile
            FileOnly(Rfileprov)
            TypY = Mid(Rfileprov, 3, 1)

            nFic = FreeFile()
            outfile = "tttemp.txt"
            If TypY = "C" Then TypY = Mid(Rfileprov, 3, 2)
            FileOpen(nFic, outfile, OpenMode.Output)
            PrintLine(nFic, NCurve, ",", TesGraph, ",", NoAxeX, ",", NoAxeY, ",", TypY)
            FileClose(nFic)

            frm01.ShowDialog()
            If frm01.Option1.Checked = True Then NoAxeX = 0
            If frm01.Option2.Checked = True Then NoAxeX = 1
            If frm01.Option8.Text = "humidit� relative" Then NoAxeY = 0
            If frm01.Option8.Text = "teneur en eau" Then NoAxeY = 1
            If frm01.Option8.Text = "teneur en ions chlorures totaux" Then NoAxeY = 2
            If frm01.Option8.Text = "teneur en ions chlorures libres" Then NoAxeY = 3
            If frm01.Option8.Text = "temp�rature" Then NoAxeY = 4
            If frm01.Option4.Checked = True Then NoAxeY = 5
            If frm01.Option8.Text = "valeur du PH" Then NoAxeY = 6
            If frm01.Option8.Text = "profondeur de carbonatation" Then NoAxeY = 7
            CoorX = CType(frm01.Text1.Text, Single)
            frm01.Close()
        Else
            frm01.ShowDialog()
            CoorX = CType(frm01.Text1.Text, Single)
            frm01.Close()
        End If
        If NoAxeY = 5 Or NoAxeY = 1 Then
            k = 2
        Else
            k = 1
        End If

        'pointage de la colonne ou de la ligne
        If NoAxeY = 7 Then
            If frm01.RadioButton1.Checked = True Then
                CoPosX = 3
            Else
                CoPosX = 4
            End If
        Else
            If NoAxeX = 0 Then 'le temps
                For i = 3 To NCol - k
                    If CoorX >= Data(1, i) And CoorX <= Data(1, i + 1) Then
                        If CoorX - Data(1, i) < Data(1, i + 1) - CoorX Then
                            CoPosX = i
                            Exit For
                        Else
                            CoPosX = i + 1
                            Exit For
                        End If
                    End If
                Next i
                If CStr(Data(1, NCol - k)) = "cumul" Then
                    If CoorX >= Data(1, NCol - k - 1) Then CoPosX = NCol - k - 1
                Else
                    If CoorX >= Data(1, NCol - k) Then CoPosX = NCol - k
                End If
                If NoAxeY = 5 Then CoPosX = NCol - k + 1
            Else 'la profondeur
                For i = 2 To NLine + 1
                    If CoorX >= Data(i, 2) And CoorX < Data(i + 1, 2) Then
                        If CoorX - Data(i, 2) < Data(i + 1, 2) - CoorX Then
                            CoPosX = i
                        Else
                            CoPosX = i + 1
                        End If
                    End If
                Next i
                If CoorX >= Data(NLine + 2, 2) Then CoPosX = NLine + 2
            End If
        End If
        cpteur = cpteur + 1

        'test la valeur maximale et minimale sur les deux axes
        If TesGraph = 0 Or TesGraph = 2 Then
            If NoAxeX = 0 Then 'valeur max et min de x
                x1Min = x2Min
                x1Max = x2Max
            End If

            If TesGraph = 0 Then
                XMI = x1Min
                XMA = x1Max
            End If

            YMI = yMin
            YMA = yMax
        End If

        outfile = "G_temp_" & NCurve & ".txt"
        nFic = FreeFile()
        FileOpen(nFic, outfile, OpenMode.Output)
        If NoAxeX = 0 Then 'le temps
            PrintLine(nFic, NLine + 1)
            PrintLine(nFic, NoAxeX)
            PrintLine(nFic, NoAxeY)
            PrintLine(nFic, CoorX)
            PrintLine(nFic, Data(2, 2))
            PrintLine(nFic, Data(2, CoPosX))
            For i = 3 To NLine + 2
                PrintLine(nFic, Data(i, 2))
                PrintLine(nFic, Data(i, CoPosX))
            Next i
        Else 'la profondeur
            PrintLine(nFic, NCol - k - 2)
            PrintLine(nFic, NoAxeX)
            PrintLine(nFic, NoAxeY)
            PrintLine(nFic, CoorX)
            PrintLine(nFic, Data(1, 3))
            PrintLine(nFic, Data(CoPosX, 3))
            For i = 4 To NCol - k
                PrintLine(nFic, Data(1, i))
                PrintLine(nFic, Data(CoPosX, i))
            Next i
        End If

        FileClose(nFic)

        TesGraph = 1

b:      'user pressed cancel error

    End Sub

    'Cherche le fichier r�sultat des points et enregistre dans un fichier texte
    Public Sub Design_point(ByRef frm02 As Chloride_transport.frmChlor, ByRef Data(,) As Object, ByRef nrofil As Short, ByRef NoAxeX As Short, ByRef NoAxeY As Short, ByRef TesPoint As Short, ByVal XMI As Single, ByRef XMA As Single, ByRef YMI As Single, ByRef YMA As Single)
        Dim Filtre As String
        Dim Index As Short
        Dim Directoire As Boolean
        Dim Titre As String
        Dim Canc As Boolean = False
        Dim OutFile As String
        Dim xMin, xMax, yMin, yMax As Single
        Dim dim1 As Short
        Dim nFic As Short
        Dim i As Short
        xMin = 9999999
        xMax = -9999999
        yMin = 9999999
        yMax = -9999999

        If TesPoint = 0 Then
            ''''''''''''''''''''''''''''''''''''''''''''
            Filtre = "txt files (P_*.txt)|P_*.txt"
            Index = 2
            Directoire = True
            Titre = "Ouvrir le fichier de profil mesur�"
            OpenDialog(OutFile, Canc, Filtre, Index, Directoire, Titre)
            If Canc = True Then GoTo b
            '''''''''''''''''''''''''''''''''''''''''''''
            LectPfData(dim1, Data, OutFile, xMin, xMax, yMin, yMax)

            If XMI > xMin Then XMI = xMin
            If XMA < xMax Then XMA = xMax
            If YMI > yMin Then YMI = yMin
            If YMA < yMax Then YMA = yMax
        End If

        OutFile = "P_temp_" & nrofil & ".txt"
        nFic = FreeFile()
        FileOpen(nFic, OutFile, OpenMode.Output, OpenAccess.Write, OpenShare.Shared)
        Print(nFic, dim1)
        For i = 1 To CInt(dim1)
            Print(nFic, Data(1, i), Data(2, i))
        Next i
        FileClose(nFic)

b:      'user pressed cancel error

    End Sub

    'Premi�re lecture du fichier de donn�es pour les points
    Public Sub LectPfData(ByRef dim1 As Short, ByRef Data As Object, ByRef OutFile As String, ByRef XMI As Single, ByRef XMA As Single, ByRef YMI As Single, ByRef YMA As Single)
        Dim nFic As Short
        Dim j As Short

        nFic = FreeFile()
        FileOpen(nFic, OutFile, OpenMode.Input, OpenAccess.Read, OpenShare.Shared)
        Input(nFic, dim1) 'nbre de coordonn�es
        ReDim Data(2, dim1)
        For j = 1 To CType(dim1, Short)
            Input(nFic, Data(1, j)) 'valeur x
            Input(nFic, Data(2, j)) 'valeur y
            'valeur des axes max et min
            If Data(2, j) < YMI Then YMI = Data(2, j)
            If Data(2, j) > YMA Then YMA = Data(2, j)
            If Data(1, j) < XMI Then XMI = Data(1, j)
            If Data(1, j) > XMA Then XMA = Data(1, j)
        Next j
        FileClose(nFic)
    End Sub

    'Ouverture de la bo�tre de dialogue des �chelles
    Public Sub Scale_info(ByRef frm02 As Chloride_transport.frmChlor, ByRef frm03 As Chloride_transport.frmScale1, ByRef Data As Object, ByRef NTCurve As Short, ByRef nTProfil As Short, ByRef TesGraph As Short, ByRef NoAxeX As Short, ByRef NoAxeY As Short, ByRef a As Single, ByRef TesPoint As Short, ByRef XMI As Single, ByRef XMA As Single, ByRef YMI As Single, ByRef YMA As Single, ByRef EctX As Single, ByRef EctY As Single)
        EctX = (XMA - XMI) / 10
        EctY = (YMA - YMI) / 10
        frm03.Text1.Text = CStr(XMI) 'min
        frm03.Text2.Text = CStr(XMA) 'max
        frm03.Text3.Text = CStr(EctX) '�cart
        frm03.Text4.Text = CStr(YMI) 'min
        frm03.Text5.Text = CStr(YMA) 'max
        frm03.Text6.Text = CStr(EctY) '�cart

        frm03.ShowDialog()
        TesGraph = 0

        EctX = CType(frm03.Text3.Text, Single) '�cart
        YMI = CType(frm03.Text4.Text, Single) 'min
        YMA = CType(frm03.Text5.Text, Single) 'max
        EctY = CType(frm03.Text6.Text, Single) '�cart
        XMA = CType(frm03.Text2.Text, Single) 'max
        XMI = CType(frm03.Text1.Text, Single) 'min
        frm03.Close()

    End Sub

    'Lecture des fichiers temporaires de points
    Public Sub readPtncs(ByRef Dim1 As Short, ByRef DataX() As Single, ByRef DataY() As Single, ByRef j As Short)
        Dim nFic As Short
        Dim OutFile As String
        Dim k As Short
        OutFile = "P_temp_" & j & ".txt"
        nFic = FreeFile()
        FileOpen(nFic, OutFile, OpenMode.Input, OpenAccess.Read, OpenShare.Shared)
        Input(nFic, Dim1) 'nbre de coordonn�es
        ReDim DataX(Dim1)
        ReDim DataY(Dim1)
        For k = 1 To Dim1
            Input(nFic, DataX(k))  'valeur x
            Input(nFic, DataY(k))  'valeur y
        Next k
        FileClose(nFic)
    End Sub

    Public Sub TakeResult(ByRef Nline As Short, ByRef DataX() As Single, ByRef DataY() As Single, ByRef CoorX As Single, ByRef kk As Short, ByRef j As Short)
        Dim nFic As Short
        Dim OutFile As String
        Dim k As Short
        Dim NoAxex As Short
        Dim NoAxey As Short

        OutFile = "G_temp_" & j & ".txt"
        nFic = FreeFile()
        FileOpen(nFic, OutFile, OpenMode.Input)
        Input(nFic, Nline)
        ReDim DataX(Nline + 1)
        ReDim DataY(Nline + 1)
        Input(nFic, NoAxex)
        Input(nFic, NoAxey)
        Input(nFic, CoorX)
        For k = 1 To Nline
            Input(nFic, DataX(k))
            Input(nFic, DataY(k))
        Next k
        FileClose(nFic)
    End Sub

    'Dessin des points et des courbes
    Public Sub GraphBMP(ByRef frm As frmChlor, ByRef nTprofil As Short, ByRef NTCurve As Short, ByRef Nline As Integer, ByRef CoorX As Single, ByRef kk As Short, ByRef NoAxex As Short, ByRef NoAxey As Short, ByRef x1 As Single, ByRef x2 As Single, ByRef h1 As Single, ByRef h2 As Single, ByRef EctX As Single, ByRef EctY As Single)
        Dim m_Image As Bitmap
        Dim m_Gr As Graphics
        Dim m_Pen As Pen
        Dim Coloras As System.Drawing.Color
        Dim Colorbs As System.Drawing.Color
        Dim DataX(1) As Single
        Dim DataY(1) As Single
        Dim i As Short
        Dim j As Integer
        Dim m_Facteur As Single

        m_Facteur = frm.Width / 17

        Dim Canc As Boolean = False
        Dim ColorDialog1 As New ColorDialog

        Dim Titel As String
        Dim Message0 As String
        Dim Message1 As String
        Dim Message2 As String
        Dim Message3 As String
        Dim Unit As String
        Dim LibelX As String
        Dim LibelY As String
        Dim msg As String

        Dim hEcart As Single

        Dim m_X1 As Single
        Dim m_X2 As Single
        Dim m_Y2 As Single
        Dim m_Y1 As Single
        Dim x As Single
        Dim y As Single

        frm.PictureBox1.Width = frm.Width
        frm.PictureBox1.Height = frm.Height
        frm.PictureBox1.Location = New Point(0, 0)
        m_Image = New Bitmap(frm.PictureBox1.Width, frm.PictureBox1.Height)
        m_Gr = Graphics.FromImage(m_Image)

        m_X1 = m_Facteur
        m_Y1 = frm.PictureBox1.Height - m_Facteur * 1.5
        m_X2 = frm.PictureBox1.Width - 0.5 * m_Facteur
        m_Y2 = m_Facteur

        Dim sh As Single = (m_Y1 - m_Y2) / (h2 - h1)
        Dim ih As Integer = CInt((h2 - h1) / EctY)

        Dim sx As Single = (m_X2 - m_X1) / (x2 - x1)
        Dim ix As Integer = CInt((x2 - x1) / EctX)

        Dim H8ft As Single = 0.2118 * CSng(m_Facteur)
        Dim H10ft As Single = 0.2648 * CSng(m_Facteur)
        Dim H12ft As Single = 0.3178 * CSng(m_Facteur)
        Dim fnt8Arial As New Font("Arial", H8ft)
        Dim fnt10Arial As New Font("Arial", H10ft)
        Dim fnt12Arial As New Font("Arial", H12ft)
        Dim pen01 As New Drawing.Pen(System.Drawing.Color.Gray, 1)
        Dim dPoint As Short = CType(0.2648 * m_Facteur, Short)
        Dim foreBrush As Brush
        Dim Drawformat0 As New StringFormat
        Dim Drawformat1 As New StringFormat
        Dim Drawformat2 As New StringFormat
        Drawformat0.Alignment = StringAlignment.Near
        Drawformat1.Alignment = StringAlignment.Far
        Drawformat2.Alignment = StringAlignment.Center

        TakeResult(Nline, DataX, DataY, CoorX, kk, 1)
        LabelAxeX(NoAxex, Unit, LibelX)
        LabelAxeY(NoAxey, Titel, hEcart, LibelY, CoorX, Unit)
        Titel = InputBox("Titre du graphique", "Editeur de graphique", Titel)
        LibelX = InputBox("L�gende axe X", "Editeur de graphique", LibelX)
        LibelY = InputBox("L�gende axe Y", "Editeur de graphique", LibelY)
        Message0 = Titel
        Message1 = LibelY
        Message2 = LibelX
        Message3 = InputBox("L�gende en-dessous du graphique", "Editeur de graphique")

        Dim drawPoint04 As New PointF(10 * m_Facteur, 0)
        m_Gr.DrawString(Message0, fnt12Arial, Brushes.Black, drawPoint04, Drawformat2)

        Dim drawPoint01 As New PointF(m_X1 - m_Facteur / 6, m_Y2 - 0.5 * m_Facteur)
        m_Gr.DrawString(Message1, fnt10Arial, Brushes.Black, drawPoint01, Drawformat0)

        Dim drawPoint03 As New PointF(m_X2, m_Y1 + 0.3 * m_Facteur)
        m_Gr.DrawString(Message2, fnt10Arial, Brushes.Black, drawPoint03, Drawformat1)

        Dim drawPoint00 As New PointF(m_Facteur * 10, frm.PictureBox1.Height - 2 * H12ft)
        m_Gr.DrawString(Message3, fnt12Arial, Brushes.Black, drawPoint00, Drawformat2)

        For i = 0 To ih
            y = m_Y1 - i * sh * EctY
            Dim drawPoint As New PointF(m_X1, y - 0.5 * H8ft)
            ' Create a Pen object
            m_Gr.DrawLine(pen01, m_X1, y, m_X2, y)
            msg = Format(h1 + CSng(i) * EctY, "0.0")
            m_Gr.DrawString(msg, fnt8Arial, Brushes.Black, drawPoint, Drawformat1)
        Next i

        For i = 0 To ix
            x = m_X1 + i * EctX * sx
            ' Create a Pen object.
            Dim drawPoint02 As New PointF(x, m_Y1)
            m_Gr.DrawLine(pen01, x, m_Y1, x, m_Y2)
            msg = i * EctX + x1
            m_Gr.DrawString(msg, fnt8Arial, Brushes.Black, drawPoint02, Drawformat2)
        Next i

        For i = 1 To NTCurve
            TakeResult(Nline, DataX, DataY, CoorX, kk, i)
            '_________________________________________________

            'Afficher la bo�te de dialogue Couleur
            With ColorDialog1
                .SolidColorOnly = True
                If .ShowDialog() = DialogResult.OK Then
                    Colorbs = ColorDialog1.Color
                    Canc = False
                Else
                    Canc = True
                End If
            End With
            If Canc = True Then GoTo b
            '_________________________________________________

            Try

                m_Pen = New Pen(Colorbs, 6)

                'dessin des r�sultats
                Dim x_old As Single
                Dim y_old As Single

                y = m_Y1 - sh * (DataY(1) - h1)
                x = m_X1 + sx * (DataX(1) - x1)
                x_old = m_X1
                y_old = y
                For j = 2 To Nline
                    y = m_Y1 - sh * (DataY(j) - h1)
                    x = m_X1 + sx * (DataX(j) - x1)
                    m_Gr.DrawLine(m_Pen, x_old, y_old, x, y)
                    x_old = x
                    y_old = y
                Next j
            Catch e As Exception
            End Try
        Next

        For i = 1 To nTprofil
            readPtncs(Nline, DataX, DataY, i)
            '_________________________________________________

            'Afficher la bo�te de dialogue Couleur
            With ColorDialog1
                .SolidColorOnly = True
                If .ShowDialog() = DialogResult.OK Then
                    Coloras = ColorDialog1.Color
                    Canc = False
                Else
                    Canc = True
                End If
            End With
            If Canc = True Then GoTo b
            '_________________________________________________

            foreBrush = New SolidBrush(Coloras)
            For j = 1 To Nline
                y = m_Y1 - sh * (DataY(j) - h1) - CShort(0.2648 * m_Facteur) / 2
                x = m_X1 + sx * (DataX(j) - x1) - CShort(0.2648 * m_Facteur) / 2
                m_Gr.FillEllipse(foreBrush, x, y, dPoint, dPoint)
            Next j
        Next

        frm.PictureBox1.Image = m_Image
b:      'user pressed cancel error
    End Sub

    'D�finition de la l�gende des axes y avec leur unit� et avec l'�cart initial dans le graphique
    Public Sub LabelAxeY(ByRef NoAxey As Short, ByRef Titel As String, ByRef hEcart As Single, ByRef LibelY As String, ByRef CoorX As Single, ByRef Unit As String)
        Select Case NoAxey
            Case 0
                Titel = "Moisture Potential Distribution at " & Format(CoorX, "0") & Unit
                hEcart = 0.1
                LibelY = "Moisture Potential [P/Ps]"
            Case 1
                Titel = "Moisture Content Distribution at " & Format(CoorX, "0") & Unit
                hEcart = 10
                LibelY = "Moisture Content [kg/m3]"
            Case 2
                Titel = "Total Chloride Ion Distribution at " & Format(CoorX, "0") & Unit
                hEcart = 0.5
                LibelY = "Total Chloride Ion Content [kg/m3]"
            Case 3
                Titel = "Free Chloride Ion Distribution at " & Format(CoorX, "0") & Unit
                hEcart = 0.1
                LibelY = "Free Chloride Ion Content [kg/m3]"
            Case 4
                Titel = "Temperature Potential Distribution at " & Format(CoorX, "0") & Unit
                hEcart = 5
                LibelY = "Temperature Potential [�C]"
            Case 5
                Titel = "Moisture Content Distribution Cumul"
                hEcart = 10
                LibelY = "Moisture Content [kg/m3]"
            Case 6
                Titel = "Evolution PH"
                hEcart = 1
                LibelY = "PH"
            Case 7
                Titel = "Evolution of Depth carbonation"
                hEcart = 1
                LibelY = "Depth carbonation [mm]"
            Case 8
                Titel = "Evolution de la s�curit� structurale (ions chlorures)"
                hEcart = 0.1
                LibelY = "Probabilit� d'initiation de la corrosion"
            Case 9
                Titel = "Evolution de la s�curit� structurale (ions chlorures)"
                hEcart = 0.1
                LibelY = "Indice de fiabilit�"
            Case 10
                Titel = "Evolution de la s�curit� structurale (carbonatation)"
                hEcart = 0.1
                LibelY = "Probabilit� d'initiation de la corrosion"
            Case 11
                Titel = "Evolution de la s�curit� structurale (carbonatation)"
                hEcart = 0.1
                LibelY = "Indice de fiabilit�"
        End Select

    End Sub

    'D�finition de la l�gende des axes x avec leur unit�
    Public Sub LabelAxeX(ByRef NoAxex As Short, ByRef Unit As String, ByRef LibelX As String)
        If NoAxex = 0 Then 'le temps
            Unit = " mm"
            LibelX = "time [days]"
        ElseIf NoAxex = 1 Then 'la profondeur
            Unit = " days"
            LibelX = "depth [mm]"
        Else
            Unit = " mm"
            LibelX = "time [years]"
        End If
    End Sub

    'Enlever les bruits du � l'approche probabiliste
    Public Sub ProbGraph(ByRef frm02 As Chloride_transport.frmChlor, ByRef frm03 As Chloride_transport.frmScale1)

        Dim Filtre As String
        Dim Index As Short
        Dim Directoire As Boolean
        Dim Titre As String
        Dim OutFile As String
        Dim PostFile As String
        Dim Canc As Boolean = False
        Dim nFic As Short

        Dim Titres(7) As String
        Dim i As Integer
        Dim j As Short
        Dim k As Short
        Dim m As Short
        Dim Data(70000, 2) As Double
        Dim Npoint(2, 2) As Double
        Dim Tijd As Double
        Dim deltaT As Double
        Dim Var As String
        Dim XMA As Double
        Dim XMI As Double
        Dim YMA As Double
        Dim YMI As Double
        Dim TestB As String
        Dim NoAxeY As Short
        Dim EctX As Single
        Dim EctY As Single
        Dim Surf As Double
        Dim NbrePoint As Integer
        Dim Msg As String
        Dim Prov As Double
        'Dim Tst As Boolean
        'Dim F1 As Worksheet
        'Dim F2 As Worksheet

        Filtre = "txt files |BCa_*.txt;BCL_*.txt;PFCa_*.txt;PFCL_*.txt|All files (*.*)|*.*"
        Index = CShort(1)
        Directoire = True
        Titre = "S�lectionner un fichier r�sultat � traiter"
        OpenDialog(OutFile, Canc, Filtre, Index, Directoire, Titre)
        If Canc = True Then GoTo b
        '''''''''''''''''''''''''''''''''''''''''''''
        nFic = FreeFile()
        FileOpen(nFic, OutFile, OpenMode.Input, OpenAccess.Read, OpenShare.Shared)
        For i = 0 To 7      'lecture du fichier
            Input(nFic, Titres(i))
        Next i
        i = 0
        YMI = 9999999999
        YMA = 0
        Do While i >= 0
            For j = 0 To 2
                Try
                    Input(nFic, Data(i, j)) 'lecture du fichier r�sultat
                    If YMI > Data(i, 2) Then YMI = Data(i, 2)
                    If YMA < Data(i, 2) Then YMA = Data(i, 2)
                Catch ex As Exception
                    Exit Do
                End Try
            Next j
            Input(nFic, Var)
            i = i + 1
        Loop
        FileClose(nFic)
        i = i - 1
        TestB = MsgBox("L'�chelle du temps en jours ?", MsgBoxStyle.YesNo, "Axe des abcisses")
        If TestB = MsgBoxResult.Yes Then
            XMI = Data(0, 1)
            XMA = Data(i, 1)
            k = 1
        Else
            XMI = Data(0, 0)
            XMA = Data(i, 0)
            k = 0
        End If
        FileOnly(OutFile)
        If Left(OutFile, 1) = "B" Then
            OutFile = Right(OutFile, Len(OutFile) - 1)
            NoAxeY = 9
            YMA = 6.0
            YMI = -6.0
        Else
            OutFile = Right(OutFile, Len(OutFile) - 2)
            NoAxeY = 8
        End If
        If Left(OutFile, 2) = "Ca" Then NoAxeY = NoAxeY + 2
        EctX = (XMA - XMI) / 10
        EctY = (YMA - YMI) / 10
c:      Scale_info(frm02, frm03, Data, 2, 0, 0, 2, NoAxeY, 0, 0, XMI, XMA, YMI, YMA, EctX, EctY)

        Dessin(frm02, 1, i, k, Msg, 2, NoAxeY, XMI, XMA, YMI, YMA, EctX, EctY, Data)
        If Msg = "" Then Msg = " "
        TestB = MsgBox("Changer l'�chelle du graphique ?", MsgBoxStyle.YesNo, "Fichier de donn�es")
        If TestB = MsgBoxResult.Yes Then GoTo c

a:      Try
            NbrePoint = InputBox("Nombre de points totaux sur le graphique", "Choix", i + 1)
        Catch ex As Exception
            GoTo b
        End Try
        Surf = 0
        Tijd = 0
        m = 0
        If i + 1 > NbrePoint Then
            ReDim Npoint(NbrePoint, 3)
            deltaT = (Data(i, 1) - Data(0, 1)) / (NbrePoint - 1)
            For j = 0 To i - 1
                Surf = 0.5 * (Data(j, 2) + Data(j + 1, 2)) * (Data(j + 1, 1) - Data(j, 1)) + Surf
                Tijd = Tijd + Data(j + 1, 1) - Data(j, 1)
                If Tijd > deltaT * m + deltaT / 2 Then
                    Npoint(m, 1) = deltaT * m
                    Npoint(m, 0) = Npoint(m, 1) / 365
                    Prov = (2 * Data(j + 1, 2) - (Data(j + 1, 2) - Data(j, 2)) * (Tijd - (deltaT * m + deltaT / 2)) / (Tijd - (Tijd - (Data(j + 1, 1) - Data(j, 1))))) * 0.5 * (Tijd - (deltaT * m + deltaT / 2))
                    Surf = Surf - Prov
                    Npoint(m, 2) = Surf / deltaT
                    If m = 0 Then Npoint(m, 2) = Surf / (deltaT / 2)
                    m = m + 1
                    Surf = Prov
                End If
            Next
            Npoint(m, 0) = Data(i, 0)
            Npoint(m, 1) = Data(i, 1)
            Npoint(m, 2) = Surf / deltaT
        End If

        Dessin(frm02, 1, NbrePoint - 1, k, Msg, 2, NoAxeY, XMI, XMA, YMI, YMA, EctX, EctY, Npoint)
        TestB = MsgBox("Changer le nombre de points ?", MsgBoxStyle.YesNo, "Graphique")
        If TestB = MsgBoxResult.Yes Then GoTo a

        TestB = MsgBox("Emp�cher un accroissement de la r�sistance ?", MsgBoxStyle.YesNo, "Graphique")
        If TestB = MsgBoxResult.Yes Then
            If NoAxeY = 9 Or NoAxeY = 11 Then
                For j = 1 To NbrePoint
                    If Npoint(j, 2) > Npoint(j - 1, 2) Then Npoint(j, 2) = Npoint(j - 1, 2)
                Next
            Else
                For j = 1 To NbrePoint
                    If Npoint(j, 2) < Npoint(j - 1, 2) Then Npoint(j, 2) = Npoint(j - 1, 2)
                Next
            End If
        End If
        Dessin(frm02, 1, NbrePoint - 1, k, Msg, 2, NoAxeY, XMI, XMA, YMI, YMA, EctX, EctY, Npoint)

        OutFile = Right(OutFile, Len(OutFile) - 3)
        OutFile = "FIN_" & OutFile        'enregistrement des donn�es
        ''''''''''''''''''''''''''''''''''
        Filtre = "txt files (FIN_*.txt)|FIN_*.txt|All files (*.*)|*.*"
        Index = CShort(1)
        Directoire = True
        Titre = "Param�tres des lois probabilistes"
        SaveDialog(OutFile, Canc, Filtre, Index, Directoire, Titre)
        If Canc = True Then GoTo d
        ''''''''''''''''''''''''''''''''''
        nFic = CShort(FreeFile())
        FileOpen(CInt(nFic), OutFile, OpenMode.Output)
        PrintLine(CInt(nFic), Titres(0))
        PrintLine(CInt(nFic), Titres(1), ",", Titres(2), ",", Titres(3), ",")
        PrintLine(CInt(nFic), Titres(5), ",", Titres(6), ",")
        For i = 0 To NbrePoint - 1
            PrintLine(CInt(nFic), Npoint(i, 0), ",", Npoint(i, 1), ",", Npoint(i, 2), ",")
        Next i
        FileClose(nFic)

d:      FileOnly(OutFile)
        OutFile = Right(OutFile, Len(OutFile) - 4)
        OutFile = Left(OutFile, Len(OutFile) - 4)
        ''''''''''''''''''''''''''''''''''
        Filtre = "txt files (GR_*.BMP)|GR_*.BMP"    'enregistrement du graph
        Index = 1
        Directoire = True
        Titre = "Enregistrer le graphique"
        SaveDialog(OutFile, Canc, Filtre, Index, Directoire, Titre)
        If Canc = True Then GoTo b

        ''''''''''''''''''''''''''''''''''
        FilePost(OutFile, PostFile)
        FileOnly(OutFile)
        OutFile = "GR_" & OutFile
        frm02.PictureBox1.Image.Save(OutFile)

b:
        frm02.PictureBox6.Visible = True
        frm02.Hide()
        frm02.Close()

    End Sub

    'Dessin des courbes probabilistiques
    Public Sub Dessin(ByRef frm As frmChlor, ByRef NTCurve As Short, ByRef Nline As Integer, ByRef CoorX As Single, ByRef Message3 As String, ByRef NoAxex As Short, ByRef NoAxey As Short, ByRef x1 As Single, ByRef x2 As Single, ByRef h1 As Single, ByRef h2 As Single, ByRef EctX As Single, ByRef EctY As Single, ByRef Data(,) As Double)
        Dim m_Image As Bitmap
        Dim m_Gr As Graphics
        Dim m_Pen As Pen
        'Dim Coloras As System.Drawing.Color
        Dim Colorbs As System.Drawing.Color
        Dim i As Short
        Dim j As Integer
        Dim m_Facteur As Single

        m_Facteur = frm.Width / 17

        Dim Canc As Boolean = False
        Dim ColorDialog1 As New ColorDialog

        Dim Titel As String
        Dim Message0 As String
        Dim Message1 As String
        Dim Message2 As String
        Dim Unit As String
        Dim LibelX As String
        Dim LibelY As String
        Dim msg As String

        Dim hEcart As Single

        Dim m_X1 As Single
        Dim m_X2 As Single
        Dim m_Y2 As Single
        Dim m_Y1 As Single
        Dim x As Single
        Dim y As Single

        frm.PictureBox1.Width = frm.Width
        frm.PictureBox1.Height = frm.Height
        frm.PictureBox1.Location = New Point(0, 0)
        m_Image = New Bitmap(frm.PictureBox1.Width, frm.PictureBox1.Height)
        m_Gr = Graphics.FromImage(m_Image)

        m_X1 = m_Facteur
        m_Y1 = frm.PictureBox1.Height - m_Facteur * 1.5
        m_X2 = frm.PictureBox1.Width - 0.5 * m_Facteur
        m_Y2 = m_Facteur

        Dim sh As Single = (m_Y1 - m_Y2) / (h2 - h1)
        Dim ih As Integer = CInt((h2 - h1) / EctY)

        Dim sx As Single = (m_X2 - m_X1) / (x2 - x1)
        Dim ix As Integer = CInt((x2 - x1) / EctX)

        Dim H8ft As Single = 0.2118 * CSng(m_Facteur)
        Dim H10ft As Single = 0.2648 * CSng(m_Facteur)
        Dim H12ft As Single = 0.3178 * CSng(m_Facteur)
        Dim fnt8Arial As New Font("Arial", H8ft)
        Dim fnt10Arial As New Font("Arial", H10ft)
        Dim fnt12Arial As New Font("Arial", H12ft)
        Dim pen01 As New Drawing.Pen(System.Drawing.Color.Gray, 1)
        Dim dPoint As Short = CType(0.2648 * m_Facteur, Short)
        Dim foreBrush As Brush
        Dim Drawformat0 As New StringFormat
        Dim Drawformat1 As New StringFormat
        Dim Drawformat2 As New StringFormat
        Drawformat0.Alignment = StringAlignment.Near
        Drawformat1.Alignment = StringAlignment.Far
        Drawformat2.Alignment = StringAlignment.Center
        If CoorX = 0 Then
            NoAxex = 2
        Else
            NoAxex = 0
        End If
        LabelAxeX(NoAxex, Unit, LibelX)
        LabelAxeY(NoAxey, Titel, hEcart, LibelY, CoorX, Unit)
        'Titel = InputBox("Titre du graphique", "Editeur de graphique", Titel)
        'LibelX = InputBox("L�gende axe X", "Editeur de graphique", LibelX)
        'LibelY = InputBox("L�gende axe Y", "Editeur de graphique", LibelY)
        Message0 = Titel
        Message1 = LibelY
        Message2 = LibelX
        If Message3 = "" Then Message3 = InputBox("L�gende en-dessous du graphique", "Editeur de graphique")

        Dim drawPoint04 As New PointF(10 * m_Facteur, 0)
        m_Gr.DrawString(Message0, fnt12Arial, Brushes.Black, drawPoint04, Drawformat2)

        Dim drawPoint01 As New PointF(m_X1 - m_Facteur / 6, m_Y2 - 0.5 * m_Facteur)
        m_Gr.DrawString(Message1, fnt10Arial, Brushes.Black, drawPoint01, Drawformat0)

        Dim drawPoint03 As New PointF(m_X2, m_Y1 + 0.3 * m_Facteur)
        m_Gr.DrawString(Message2, fnt10Arial, Brushes.Black, drawPoint03, Drawformat1)

        Dim drawPoint00 As New PointF(m_Facteur * 10, frm.PictureBox1.Height - 2 * H12ft)
        m_Gr.DrawString(Message3, fnt12Arial, Brushes.Black, drawPoint00, Drawformat2)

        For i = 0 To ih
            y = m_Y1 - i * sh * EctY
            Dim drawPoint As New PointF(m_X1, y - 0.5 * H8ft)
            ' Create a Pen object
            m_Gr.DrawLine(pen01, m_X1, y, m_X2, y)
            msg = Format(h1 + CSng(i) * EctY, "0.0")
            m_Gr.DrawString(msg, fnt8Arial, Brushes.Black, drawPoint, Drawformat1)
        Next i

        For i = 0 To ix
            x = m_X1 + i * EctX * sx
            ' Create a Pen object.
            Dim drawPoint02 As New PointF(x, m_Y1)
            m_Gr.DrawLine(pen01, x, m_Y1, x, m_Y2)
            msg = i * EctX + x1
            m_Gr.DrawString(msg, fnt8Arial, Brushes.Black, drawPoint02, Drawformat2)
        Next i

        For i = 1 To NTCurve
            ''_________________________________________________

            ''Afficher la bo�te de dialogue Couleur
            'With ColorDialog1
            '    .SolidColorOnly = True
            '    If .ShowDialog() = DialogResult.OK Then
            '        Colorbs = ColorDialog1.Color
            '        Canc = False
            '    Else
            '        Canc = True
            '    End If
            'End With
            'If Canc = True Then GoTo b
            ''_________________________________________________

            Try
                Colorbs = System.Drawing.Color.Black
                m_Pen = New Pen(Colorbs, 6)

                'dessin des r�sultats
                Dim x_old As Single
                Dim y_old As Single

                y = m_Y1 - sh * (Data(0, 2) - h1)
                x = m_X1 + sx * (Data(0, CoorX) - x1)
                x_old = m_X1
                y_old = y
                For j = 1 To Nline
                    y = m_Y1 - sh * (Data(j, 2) - h1)
                    x = m_X1 + sx * (Data(j, CoorX) - x1)
                    m_Gr.DrawLine(m_Pen, x_old, y_old, x, y)
                    x_old = x
                    y_old = y
                Next j
            Catch e As Exception
            End Try
        Next

        frm.PictureBox1.Image = m_Image
b:      'user pressed cancel error
    End Sub

End Module
