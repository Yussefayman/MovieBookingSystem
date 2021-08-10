package moviebookingsystem;

import database.mysql;
import java.io.IOException;
import java.net.URL;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ResourceBundle;
import javafx.event.ActionEvent;
import javafx.fxml.FXML;
import javafx.fxml.Initializable;
import javafx.scene.control.Button;
import javafx.scene.control.Label;


public class Poster_infoController implements Initializable {

    @FXML Label moviename;
    @FXML Label moviegenre;
    @FXML Label moviedesc;
    @FXML Label movierele_date;
    @FXML Button button ;
    
    @FXML void showDesc(ActionEvent e ) throws SQLException{
    
        int movieId= Integer.parseInt( moviename.getText());
       
        ResultSet rs = mysql.RunOnQuery("SELECT * FROM movie where movie_id = '"+movieId+"' ");
        while (rs.next()){
            String name = rs.getString("movie_name");
           
            String genre = rs.getString("movie_genre");
            String movieD = rs.getString("movie_desc");
            String releD = rs.getString("movie_release_date");
            moviegenre.setText(genre);
            moviename.setVisible(true);
            moviename.setText(name);
            moviedesc.setText(movieD);
            movierele_date.setText(releD);
            
            break;
        }
    }
    
    @FXML
    public void book_film(ActionEvent event) throws IOException{
        
        Controls control = new Controls();
         control.switchScenes("date_time_submit");
    }
   public void setMOVIEID(int Value){
       this.moviename.setText(Integer.toString(Value));
   }
    @Override
    public void initialize(URL url, ResourceBundle rb) {
      
    }    
    
}
